#include "configuration.hpp"
#include "geometry_validation.hpp"
#include <algorithm>
#include <charconv>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <type_traits>

namespace sdbd::cli {
namespace {
// Remove surrounding whitespace without changing quoted contents.
std::string trim(std::string s) {
    const auto first=s.find_first_not_of(" \t\r\n"),last=s.find_last_not_of(" \t\r\n");
    return first==std::string::npos ? "" : s.substr(first,last-first+1);
}
// Decode an optional quoted value, rejecting trailing text and empty strings.
std::string decode(std::string s) {
    s=trim(std::move(s));
    if (s.empty()) throw std::invalid_argument("Missing value");
    if (s.front()!='\"') return s;
    std::istringstream in(s);
    std::string value;
    if (!(in>>std::quoted(value))) throw std::invalid_argument("Unterminated quoted value");
    in>>std::ws;
    if (!in.eof() || value.empty()) throw std::invalid_argument("Invalid quoted value");
    return value;
}
// Parse a complete finite numeric value, preserving double precision for phi.
template<class T> T number(const std::string& text) {
    T value{};
    const auto r=std::from_chars(text.data(),text.data()+text.size(),value);
    if (r.ec!=std::errc{} || r.ptr!=text.data()+text.size()) throw std::invalid_argument("Invalid number: "+text);
    if constexpr (std::is_floating_point_v<T>)
        if (!std::isfinite(value)) throw std::invalid_argument("Nonfinite number");
    return value;
}
// Parse an explicit boolean, without accepting partial or ambiguous spellings.
bool boolean(const std::string& text) {
    if (text=="true") return true;
    if (text=="false") return false;
    throw std::invalid_argument("Expected true or false");
}
// Return exactly three finite, positive Cartesian lengths from a box field.
Vec3 dimensions(const std::string& text) {
    std::istringstream in(text);
    in.imbue(std::locale::classic());
    long double x,y,z;
    std::string extra;
    if (!(in>>x>>y>>z) || in>>extra || !std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z) || x<=0 || y<=0 || z<=0)
        throw std::invalid_argument("box requires three positive finite lengths: x y height");
    return Vec3(x,y,z);
}
// Return a normalized absolute path relative to the declaring configuration.
std::filesystem::path path(const Value& v) {
    return std::filesystem::absolute(v.base/v.text).lexically_normal();
}
// Return metres per declared physical length unit; legacy is deliberately unconverted.
long double unit(const std::string& name) {
    if (name=="mm") return .001L;
    if (name=="cm") return .01L;
    if (name=="m") return 1;
    if (name=="legacy") return 1;
    throw std::invalid_argument("units must be mm, cm, m or legacy");
}
// Reject a nonpositive/nonfinite physical dimension or scale.
void positive(long double value,const char* name) {
    if (!std::isfinite(value) || value<=0) throw std::invalid_argument(std::string(name)+" must be positive and finite");
}
}
Fields read_configuration(const std::filesystem::path& filename) {
    std::ifstream input(filename);
    if (!input) throw std::runtime_error("Cannot open configuration: "+filename.string());
    Fields fields;
    std::string line;
    std::size_t lineno=0;
    const auto base=std::filesystem::absolute(filename).parent_path();
    while (std::getline(input,line)) {
        ++lineno;
        try {
            bool quoted=false,escaped=false;
            for (std::size_t i=0;i<line.size();++i) {
                const char c=line[i];
                if (escaped) {
                    escaped=false;
                    continue;
                }
                if (quoted && c=='\\') {
                    escaped=true;
                    continue;
                }
                if (c=='\"') quoted=!quoted;
                if (c=='#' && !quoted) {
                    line.resize(i);
                    break;
                }
            }
            if (quoted) throw std::invalid_argument("Unterminated quoted value");
            line=trim(std::move(line));
            if (line.empty()) continue;
            const auto equal=line.find('=');
            if (equal==std::string::npos) throw std::invalid_argument("Expected key = value");
            const auto key=trim(line.substr(0,equal));
            if (!fields.emplace(key,Value{decode(line.substr(equal+1)),base,filename.string()+":"+std::to_string(lineno)}).second)
                throw std::invalid_argument("Duplicate key: "+key);
        } catch (const std::exception& e) {
            throw std::invalid_argument(filename.string()+":"+std::to_string(lineno)+": "+e.what());
        }
    }
    if (input.bad()) throw std::runtime_error("Cannot read configuration");
    return fields;
}
Configuration resolve(const Fields& fields,bool generation,bool require_output) {
    geometry_detail::ValidationEnvironment environment;
    Configuration c;
    if (generation) c.generation.emplace();
    static const std::set<std::string> common={"format","units","output","quiet"};
    static const std::set<std::string> run={
        "particle","particle_checksum","box","count","seed","rng","resolution",
        "size","length","radius","spheres","phi","redrop","injection",
        "legacy_output","stl_poses"
    };
    static const std::set<std::string> prepare={
        "stl","source_units","scale","voxel_resolution","maximum_spheres","accuracy",
        "minimum_radius_voxels","minimum_center_distance_relative","padding",
        "search_window","persistence","radius_offset_voxels","confine_mesh",
        "prune_isolated_spheres","density","maximum_voxels"
    };
    if (generation && (!fields.contains("stl") || !fields.contains("units")))
        throw std::invalid_argument("generate requires stl and declared units");
    if (require_output && !fields.contains("output")) throw std::invalid_argument("Set output in the configuration or use --output");
    for (const auto& [key,v]:fields) {
        try {
            if (!common.contains(key) && !(generation ? prepare : run).contains(key))
                throw std::invalid_argument("Unknown or misplaced key: "+key+(key=="stl" ? "; prepare STL with deposit generate first" : ""));
            const auto& s=v.text;
            auto& p=c.simulation;
            if (key=="format") {
                if (s!="SALTED_2") throw std::invalid_argument("Expected format = SALTED_2; see examples");
            }
            else if (key=="units") c.units=s;
            else if (key=="output") {
                c.output=path(v);
                if (c.output.filename().empty() && c.output!=c.output.root_path()) c.output=c.output.parent_path();
            }
            else if (key=="quiet") c.quiet=boolean(s);
            else if (key=="particle") c.particle=s;
            else if (key=="particle_checksum") c.particle_checksum=s;
            else if (key=="box") c.box=dimensions(s);
            else if (key=="count") p.count=number<int>(s);
            else if (key=="seed") p.seed=number<int>(s);
            else if (key=="rng") p.rng=s;
            else if (key=="resolution") p.resolution=number<int>(s);
            else if (key=="size") c.size=number<long double>(s);
            else if (key=="length") c.length=number<long double>(s);
            else if (key=="radius") c.radius=number<long double>(s);
            else if (key=="spheres") c.spheres=number<std::size_t>(s);
            else if (key=="phi") p.phi=number<double>(s);
            else if (key=="redrop") p.redrop=boolean(s);
            else if (key=="injection") {
                if (s=="uniform") p.injection=InjectionMode::uniform;
                else if (s=="circle") p.injection=InjectionMode::circle;
                else throw std::invalid_argument("injection must be uniform or circle");
            }
            else if (key=="legacy_output") c.legacy_output=boolean(s);
            else if (key=="stl_poses") c.stl_poses=boolean(s);
            else {
                auto& g=*c.generation;
                auto& m=g.options;
                if (key=="stl") g.stl=path(v);
                else if (key=="source_units") m.source_units=s;
                else if (key=="scale") g.scale=number<long double>(s);
                else if (key=="voxel_resolution") m.voxel_resolution=number<int>(s);
                else if (key=="maximum_spheres") m.maximum_spheres=number<int>(s);
                else if (key=="accuracy") m.reconstruction_accuracy=number<float>(s);
                else if (key=="minimum_radius_voxels") m.minimum_radius_voxels=number<int>(s);
                else if (key=="minimum_center_distance_relative") m.minimum_center_distance_relative=number<float>(s);
                else if (key=="padding") m.padding=number<int>(s);
                else if (key=="search_window") m.search_window=number<int>(s);
                else if (key=="persistence") m.persistence=number<int>(s);
                else if (key=="radius_offset_voxels") m.radius_offset_voxels=number<float>(s);
                else if (key=="confine_mesh") m.confine_mesh=boolean(s);
                else if (key=="prune_isolated_spheres") m.prune_isolated_spheres=boolean(s);
                else if (key=="density") m.density=number<long double>(s);
                else if (key=="maximum_voxels") m.maximum_voxels=number<std::size_t>(s);
            }
        } catch (const std::exception& e) {
            throw std::invalid_argument(v.location+": "+e.what());
        }
    }
    unit(c.units);
    if (generation) {
        auto& g=*c.generation;
        positive(g.scale,"scale");
        if (!fields.contains("source_units")) g.options.source_units=c.units;
        if ((c.units=="legacy")!=(g.options.source_units=="legacy")) throw std::invalid_argument("Cannot infer conversion to/from legacy units");
        g.options.source_scale=g.scale*unit(g.options.source_units)/unit(c.units);
        g.options.target_units=c.units;
    } else {
        auto& p=c.simulation;
        if (p.count<=0) throw std::invalid_argument("count must be positive (number of rigid particles)");
        if (p.resolution<=0) throw std::invalid_argument("resolution must be positive (built-in discretization)");
        if (p.phi<=0) throw std::invalid_argument("phi must be positive (rotation angle in radians)");
        if (c.particle=="cube" && !fields.contains("resolution")) p.resolution=3;
        if (c.particle=="cube" && p.resolution<2) throw std::invalid_argument("Cube resolution must be at least 2");
        p.shape=c.particle=="spiky" ? 1 : c.particle=="cube" ? 2 : 0;
        const bool builtin=c.particle=="sphere" || c.particle=="spiky" || c.particle=="cube";
        if (builtin) {
            const long double nominal=60.L*(p.shape==0 ? 1 : p.shape==1 ? 2.L*p.resolution-1 : p.resolution);
            if (!fields.contains("size")) c.size=nominal;
            positive(c.size,"size");
            p.particle_scale=c.size/nominal;
        } else if (c.particle!="stick") {
            c.particle_file=path(fields.at("particle"));
            if (fields.contains("size") || fields.contains("resolution")) throw std::invalid_argument("File particles retain their prepared dimensions; remove size and resolution from the run configuration");
        }
        if (c.particle!="stick" && (fields.contains("length") || fields.contains("radius") || fields.contains("spheres")))
            throw std::invalid_argument("length, radius and spheres apply only to stick");
        if (c.particle=="stick" && (fields.contains("size") || fields.contains("resolution")))
            throw std::invalid_argument("Stick dimensions use length, radius and spheres");
        p.rng=resolve_rng_name(p.rng);
    }
    geometry_detail::arithmetic_range();
    return c;
}
void configure_domain(Configuration& c,const std::optional<ParticleGeometry>& geometry) {
    geometry_detail::ValidationEnvironment environment;
    if (geometry && geometry->units()!="legacy" && geometry->units()!=c.units)
        throw std::invalid_argument("Prepared particle units differ from simulation units; convert geometry explicitly");
    const auto x=c.box.x(),y=c.box.y(),h=c.box.z();
    if (std::floor(x)!=x || std::floor(y)!=y || x>std::numeric_limits<int>::max() || y>std::numeric_limits<int>::max())
        throw std::invalid_argument("Horizontal box lengths must be positive integer values in the declared units; no rounding is applied");
    const int cell=std::gcd(120,std::gcd(static_cast<int>(x),static_cast<int>(y)));
    const long double extent=geometry ? geometry->bounding_radius() : 60.L*c.simulation.resolution*c.simulation.particle_scale;
    const auto max_radius=geometry ? geometry->maximum_radius() : 30*c.simulation.particle_scale;
    const auto margin=std::max(960.L,2*extent+2*std::max(34.L,1.1L*max_radius)+cell);
    const auto nz=std::ceil((h+margin)/cell);
    if (!std::isfinite(nz) || nz>std::numeric_limits<int>::max() || (x/cell)*(y/cell)*nz>20'000'000)
        throw std::invalid_argument("Box requires too many search cells; use commensurate horizontal dimensions (multiples of 120 give the coarsest grid)");
    c.simulation.grid={static_cast<int>(x/cell),static_cast<int>(y/cell),static_cast<int>(nz),
        static_cast<long double>(cell),static_cast<long double>(cell),static_cast<long double>(cell),0,0,h-nz*cell};
    c.simulation.height_limit=h;
    validate_simulation(c.simulation,geometry);
    geometry_detail::arithmetic_range();
}
void write_configuration(const std::filesystem::path& filename,const Configuration& c) {
    std::ofstream out(filename);
    out.imbue(std::locale::classic());
    out<<std::setprecision(std::numeric_limits<long double>::max_digits10)<<std::boolalpha;
    out<<"format = SALTED_2\nunits = "<<c.units<<"\nparticle = ";
    auto particle=c.particle_file ? c.particle_file->lexically_relative(std::filesystem::absolute(filename).parent_path()).generic_string() : c.particle;
    if (particle.empty() && c.particle_file) particle=c.particle_file->generic_string();
    if (!c.particle_checksum.empty()) out<<std::quoted(particle)<<"\nparticle_checksum = "<<std::quoted(c.particle_checksum);
    else out<<std::quoted(particle);
    out<<"\nbox = "<<c.box<<"\ncount = "<<c.simulation.count
       <<"\nseed = "<<c.simulation.seed<<"\nrng = "<<c.simulation.rng;
    if (!c.particle_file && c.particle!="stick") out<<"\nresolution = "<<c.simulation.resolution<<"\nsize = "<<c.size;
    if (c.particle=="stick") out<<"\nlength = "<<c.length<<"\nradius = "<<c.radius<<"\nspheres = "<<c.spheres;
    out<<"\ninjection = "<<(c.simulation.injection==InjectionMode::circle ? "circle" : "uniform");
    out<<"\nphi = "<<c.simulation.phi<<"\nredrop = "<<c.simulation.redrop<<"\nquiet = "<<c.quiet
       <<"\nlegacy_output = "<<c.legacy_output<<"\nstl_poses = "<<c.stl_poses<<"\noutput = "<<std::quoted(c.output.generic_string())<<'\n';
    out.close();
    if (!out) throw std::runtime_error("Cannot write resolved configuration: "+filename.string());
}
}
