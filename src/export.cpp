#include "sdbd/export.hpp"
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>

namespace sdbd {
namespace {
// Return a classic-locale checked output stream preserving long-double text digits.
std::ofstream output(const std::filesystem::path& path) {
    std::ofstream out;out.exceptions(std::ios::failbit|std::ios::badbit);
    try { out.open(path); }
    catch (const std::ios_base::failure& error) { throw std::runtime_error("Cannot open export "+path.string()+": "+error.what()); }
    out.imbue(std::locale::classic());
    out<<std::setprecision(std::numeric_limits<long double>::max_digits10);
    return out;
}
// Reject nonfinite Cartesian data before it reaches a VTK point record.
void finite(Vec3 p) {
    if (!std::isfinite(p.x()) || !std::isfinite(p.y()) || !std::isfinite(p.z()))
        throw std::invalid_argument("Nonfinite export coordinate");
}
// Write the legacy PolyData preamble; VTK double is the portable reader type.
void header(std::ostream& out,std::size_t count,const char* title) {
    out<<"# vtk DataFile Version 3.0\n"<<title<<"\nASCII\nDATASET POLYDATA\nPOINTS "<<count<<" double\n";
}
// Write one vertex per point, followed by the point-attribute count.
void vertices(std::ostream& out,std::size_t count) {
    if (count>std::numeric_limits<int>::max()/2) throw std::invalid_argument("Legacy VTK vertex count exceeds supported index range");
    out<<"VERTICES "<<count<<' '<<2*count<<'\n';
    for (std::size_t i=0;i<count;++i) out<<"1 "<<i<<'\n';
    out<<"POINT_DATA "<<count<<'\n';
}
}
void write_deposition_vtk(const std::filesystem::path& path,std::span<const Vec3> centers,
    std::span<const long double> radii,std::span<const std::size_t> ids) {
    if (centers.size()!=radii.size() || centers.size()!=ids.size()) throw std::invalid_argument("Misaligned sphere export arrays");
    for (std::size_t i=0;i<centers.size();++i) {
        finite(centers[i]);
        if (!std::isfinite(radii[i]) || radii[i]<=0 || ids[i]>std::numeric_limits<int>::max() || (i && ids[i]<ids[i-1]))
            throw std::invalid_argument("Invalid radius or ungrouped particle identity in export");
    }
    auto out=output(path);header(out,centers.size(),"SALTED deposited spheres; solver length units");
    for (const auto& p:centers) out<<p<<'\n';
    vertices(out,centers.size());
    out<<"FIELD PointAttributes 3\nradius 1 "<<centers.size()<<" double\n";
    for (auto radius:radii) out<<radius<<'\n';
    out<<"particle_id 1 "<<centers.size()<<" int\n";
    for (auto id:ids) out<<id<<'\n';
    out<<"sphere_id 1 "<<centers.size()<<" int\n";
    std::size_t local=0;
    for (std::size_t i=0;i<ids.size();++i) { if (i==0 || ids[i]!=ids[i-1]) local=0;out<<local++<<'\n'; }
    out.close();
}
std::string file_checksum(const std::filesystem::path& path) {
    std::ifstream source(path,std::ios::binary);
    if (!source) throw std::invalid_argument("Cannot open source file: "+path.string());
    std::uint64_t hash=14695981039346656037ULL;char byte;
    while (source.get(byte)) { hash^=static_cast<unsigned char>(byte);hash*=1099511628211ULL; }
    if (!source.eof()) throw std::runtime_error("Cannot checksum source file: "+path.string());
    std::ostringstream checksum;checksum<<"fnv1a64:"<<std::hex<<std::setw(16)<<std::setfill('0')<<hash;
    return checksum.str();
}
void validate_stl_source(const ParticleGeometry& geometry) {
    if (!geometry.stl_source() || !geometry.physical_properties())
        throw std::invalid_argument("STL pose export requires retained original STL/frame provenance; regenerate older particle snapshots from STL");
    const auto& p=*geometry.physical_properties();
    if (!p.backend_version.starts_with("GEMSS ") || p.parameters.find("compute_physics=2")==std::string::npos)
        throw std::invalid_argument("STL pose export requires GEMSS target-derived physical properties");
    if (file_checksum(geometry.stl_source()->path)!=p.checksum) throw std::invalid_argument("Original STL checksum differs from particle provenance");
}
void write_stl_poses(const std::filesystem::path& path,std::span<const Vec3> centers,
    std::span<const Quaternion> orientations,const ParticleGeometry& geometry) {
    validate_stl_source(geometry);
    if (centers.size()!=orientations.size()) throw std::invalid_argument("Misaligned particle pose arrays");
    const auto& frame=*geometry.stl_source();const auto& physics=*geometry.physical_properties();
    const auto source_center=frame.orientation.inverse().apply(geometry.center_of_mass()-frame.origin);
    auto out=output(path);header(out,centers.size(),"SALTED original STL origins; WXYZ quaternions; solver length units");
    for (std::size_t i=0;i<centers.size();++i) {
        const auto origin=centers[i]-orientations[i].normalized().apply(source_center);finite(origin);out<<origin<<'\n';
    }
    vertices(out,centers.size());
    out<<"FIELD PointAttributes 3\nparticle_id 1 "<<centers.size()<<" int\n";
    for (std::size_t i=0;i<centers.size();++i) out<<i<<'\n';
    // FIELD preserves all arrays even with a reader's ReadAllScalars disabled.
    // VECTORS cannot represent a quaternion: that record has exactly 3 components.
    out<<"orientation 4 "<<centers.size()<<" double\n";
    for (const auto& rotation:orientations) { const auto q=rotation.normalized();out<<q.w<<' '<<q.x<<' '<<q.y<<' '<<q.z<<'\n'; }
    out<<"scale 1 "<<centers.size()<<" double\n";
    for (std::size_t i=0;i<centers.size();++i) out<<physics.source_scale<<'\n';
    out.close();

}
}
