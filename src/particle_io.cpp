#include "sdbd/particle.hpp"
#include "geometry_validation.hpp"
#include <charconv>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>

namespace sdbd {
namespace {
// Read the next field name from input; throw unless it equals expected.
void label(std::istream& input, const char* expected) {
    std::string actual;
    if (!(input>>actual) || actual!=expected)
        throw std::runtime_error(std::string("Expected particle field: ")+expected);
}
// Read a finite long-double particle value from input or throw on malformed/nonfinite data.
long double number(std::istream& input) {
    long double value;
    if (!(input>>value) || !std::isfinite(value)) throw std::runtime_error("Invalid particle number");
    return value;
}
// Read Cartesian x,y,z in file order and return the coordinate in declared particle units.
Vec3 point(std::istream& input) {
    // Separate statements preserve file order regardless of argument evaluation order.
    const auto x=number(input), y=number(input), z=number(input);
    return Vec3(x,y,z);
}
}
ParticleGeometry read_particle(std::istream& input,const std::filesystem::path& base) {
    geometry_detail::ValidationEnvironment environment;
    if (!input) throw std::runtime_error("Cannot read failed particle stream");
    // Use a local classic-locale reader without changing the caller's stream state.
    std::ostringstream text;
    text<<input.rdbuf();
    if (input.bad() || text.bad()) throw std::runtime_error("Cannot read particle stream");
    std::istringstream file(text.str());
    file.imbue(std::locale::classic());
    label(file,"SDBD_MULTISPHERE");
    std::string version; file>>version;
    if (version!="1" && version!="2" && version!="3") throw std::runtime_error("Unsupported particle version");
    label(file,"units"); std::string units;file>>units;
    label(file,"scale"); const auto scale=number(file);
    if (scale<=0) throw std::runtime_error("Particle scale must be positive");
    label(file,"center");
    std::string mode;
    if (!(file>>mode)) throw std::runtime_error("Missing particle center mode");
    std::optional<Vec3> center;
    if (mode=="explicit") center=point(file);
    else if (mode!="auto") throw std::runtime_error("Unknown particle center mode");
    std::optional<PhysicalProperties> properties;
    if (version!="1") {
        properties.emplace(); auto& p=*properties;
        label(file,"target_volume"); p.volume=number(file);
        label(file,"density"); p.density=number(file);
        label(file,"inertia_per_density"); for (auto& value:p.inertia_per_density) value=number(file);
        auto text_field=[&](const char* key,std::string& value) {
            label(file,key); if (!(file>>std::quoted(value))) throw std::runtime_error("Invalid provenance string");
        };
        text_field("source",p.source); text_field("checksum",p.checksum);
        text_field("backend",p.backend_version); text_field("parameters",p.parameters);
        text_field("source_units",p.source_units);
        label(file,"source_scale"); p.source_scale=number(file);
        label(file,"voxel_size"); p.voxel_size=number(file);
    }
    std::optional<StlSource> source;
    if (version=="3") {
        source.emplace();std::string path;
        label(file,"source_stl");if (!(file>>std::quoted(path)) || path.empty()) throw std::runtime_error("Invalid source STL path");
        source->path=std::filesystem::absolute(base/path).lexically_normal();
        label(file,"source_origin");source->origin=point(file);
        label(file,"source_orientation_wxyz");
        auto& q=source->orientation;q.w=number(file);q.x=number(file);q.y=number(file);q.z=number(file);
    }
    label(file,"spheres");
    std::string token;
    if (!(file>>token)) throw std::runtime_error("Missing sphere count");
    std::size_t count{};
    const auto parsed=std::from_chars(token.data(),token.data()+token.size(),count);
    if (parsed.ec!=std::errc{} || parsed.ptr!=token.data()+token.size() || count==0 || count>ParticleGeometry::maximum_spheres)
        throw std::runtime_error("Invalid sphere count");
    std::vector<Sphere> spheres;
    for (std::size_t i=0; i<count; ++i) {
        const auto p=point(file);
        const auto radius=number(file);
        spheres.push_back({p,radius});
    }
    if (file>>token) throw std::runtime_error("Unexpected trailing particle data");
    geometry_detail::arithmetic_range();
    auto geometry=ParticleGeometry(std::move(spheres),center,std::move(properties),std::move(source),std::move(units));
    return scale==1 ? geometry : geometry.scaled(scale);
}

void write_particle(std::ostream& output, const ParticleGeometry& particle,const std::filesystem::path& base) {
    std::ostringstream file;
    file.imbue(std::locale::classic());
    file<<std::setprecision(std::numeric_limits<long double>::max_digits10);
    file<<"SDBD_MULTISPHERE "<<(particle.stl_source() ? 3 : particle.physical_properties() ? 2 : 1)<<"\nunits "<<particle.units()<<"\nscale 1\ncenter ";
    if (particle.has_explicit_center()) file<<"explicit "<<particle.center_of_mass()<<'\n';
    else file<<"auto\n";
    if (particle.physical_properties()) {
        const auto& p=*particle.physical_properties();
        file<<"target_volume "<<p.volume<<"\ndensity "<<p.density<<"\ninertia_per_density";
        for (auto value:p.inertia_per_density) file<<' '<<value;
        file<<"\nsource "<<std::quoted(p.source)<<"\nchecksum "<<std::quoted(p.checksum)
            <<"\nbackend "<<std::quoted(p.backend_version)<<"\nparameters "<<std::quoted(p.parameters)
            <<"\nsource_units "<<std::quoted(p.source_units)<<"\nsource_scale "<<p.source_scale
            <<"\nvoxel_size "<<p.voxel_size<<'\n';
    }
    if (particle.stl_source()) {
        const auto& source=*particle.stl_source();auto path=source.path;
        if (!base.empty()) { const auto relative=path.lexically_relative(base);if (!relative.empty()) path=relative; }
        const auto& q=source.orientation;
        file<<"source_stl "<<std::quoted(path.generic_string())<<"\nsource_origin "<<source.origin
            <<"\nsource_orientation_wxyz "<<q.w<<' '<<q.x<<' '<<q.y<<' '<<q.z<<'\n';
    }
    file<<"spheres "<<particle.spheres().size()<<'\n';
    for (const auto& sphere : particle.spheres()) file<<sphere.center<<' '<<sphere.radius<<'\n';
    output<<file.str();
    if (!output) throw std::runtime_error("Cannot write particle stream");
}

ParticleGeometry load_particle(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot open particle file: "+path.string());
    try { return read_particle(input,std::filesystem::absolute(path).parent_path()); }
    catch (const std::exception& error) { throw std::runtime_error(path.string()+": "+error.what()); }
}
void save_particle(const std::filesystem::path& path, const ParticleGeometry& particle) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("Cannot open particle output: "+path.string());
    write_particle(output,particle,std::filesystem::absolute(path).parent_path());
    output.close();
    if (!output) throw std::runtime_error("Cannot close particle output: "+path.string());
}
}
