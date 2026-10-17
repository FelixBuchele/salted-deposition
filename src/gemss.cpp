#include "sdbd/gemss.hpp"
#include "geometry_validation.hpp"
#include <bit>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <map>
#include <queue>
#include <sstream>
#include <Eigen/Dense>
static_assert(EIGEN_VERSION_AT_LEAST(3,4,0), "GEMSS bundled libigl requires Eigen >= 3.4");
#include "GEMSS-interface.h"

namespace sdbd {
namespace {
// Reject a failed mesh/backend precondition with invalid_argument and its technical message.
void require(bool condition,const char* message) { if (!condition) throw std::invalid_argument(message); }
// Validate binary STL size/finite vertices on a little-endian host and return its FNV-1a checksum; throw on malformed input.
std::string preflight(const std::filesystem::path& path) {
    require(std::endian::native==std::endian::little,"GEMSS binary STL loader requires little-endian host");
    std::ifstream file(path,std::ios::binary);
    require(bool(file),"Cannot open binary STL");
    std::array<char,80> header{}; std::uint32_t count{};
    file.read(header.data(),80); file.read(reinterpret_cast<char*>(&count),4);
    require(bool(file) && count>0 && count<=1'000'000,"Invalid or empty binary STL triangle count");
    require(std::filesystem::file_size(path)==84ULL+50ULL*count,"Binary STL length does not match triangle count");
    for (std::uint32_t i=0;i<count;++i) {
        std::array<float,12> triangle{}; std::uint16_t attribute;
        file.read(reinterpret_cast<char*>(triangle.data()),48); file.read(reinterpret_cast<char*>(&attribute),2);
        require(bool(file),"Truncated binary STL");
        for (int j=3;j<12;++j) require(std::isfinite(triangle[j]),"Nonfinite STL vertex");
    }
    file.clear(); file.seekg(0);
    std::uint64_t hash=14695981039346656037ULL; char byte;
    while (file.get(byte)) { hash^=static_cast<unsigned char>(byte); hash*=1099511628211ULL; }
    require(file.eof(),"Cannot checksum STL");
    std::ostringstream out; out<<"fnv1a64:"<<std::hex<<std::setw(16)<<std::setfill('0')<<hash;
    return out.str();
}
// Reject degenerate triangles, open/nonmanifold edges, inconsistent winding, disconnected surfaces or zero signed volume; self-intersections remain unchecked.
void validate_mesh(const GEMSS::STLMesh& mesh) {
    require(!mesh.is_empty() && mesh.vertices.allFinite(),"Empty or nonfinite STL mesh");
    // Incident face count/orientations and adjacency for a welded mesh edge.
    struct Edge { int count{},orientation{}; std::vector<int> faces; };
    std::map<std::pair<int,int>,Edge> edges;
    std::vector<std::vector<int>> adjacent(mesh.triangles.rows());
    const Eigen::Vector3d origin=mesh.vertices.row(0).cast<double>();
    double volume6=0;
    for (int i=0;i<mesh.triangles.rows();++i) {
        const Eigen::Vector3d a=mesh.vertices.row(mesh.triangles(i,0)).cast<double>();
        const Eigen::Vector3d b=mesh.vertices.row(mesh.triangles(i,1)).cast<double>();
        const Eigen::Vector3d c=mesh.vertices.row(mesh.triangles(i,2)).cast<double>();
        require((b-a).cross(c-a).squaredNorm()>0,"Degenerate STL triangle (possibly after scaling)");
        volume6+=(a-origin).dot((b-origin).cross(c-origin));
        for (int k=0;k<3;++k) {
            int x=mesh.triangles(i,k),y=mesh.triangles(i,(k+1)%3);
            auto& edge=edges[{std::min(x,y),std::max(x,y)}];
            ++edge.count; edge.orientation+=x<y ? 1 : -1; edge.faces.push_back(i);
        }
    }
    for (const auto& [key,edge]:edges) {
        (void)key;
        require(edge.count==2 && edge.orientation==0,"STL must be closed, manifold and consistently oriented");
        adjacent[edge.faces[0]].push_back(edge.faces[1]); adjacent[edge.faces[1]].push_back(edge.faces[0]);
    }
    std::vector<bool> seen(adjacent.size()); std::queue<int> queue; queue.push(0); seen[0]=true;
    std::size_t visited=0;
    while (!queue.empty()) { int i=queue.front();queue.pop();++visited;
        for (auto j:adjacent[i]) if (!seen[j]) { seen[j]=true;queue.push(j); }
    }
    require(visited==adjacent.size(),"Disconnected or nested STL shells are not supported by this adapter");
    require(std::isfinite(volume6) && std::abs(volume6)>0,"STL has no meaningful enclosed volume");
}
}
namespace {
// Load and validate scaled mesh/options without voxelization; return the checked mesh.
GEMSS::STLMesh checked_mesh(const std::filesystem::path& path,const GemssOptions& o,std::string* checksum=nullptr) {
    geometry_detail::ValidationEnvironment environment;
    require(o.voxel_resolution>=4 && o.voxel_resolution<=2048 && o.padding>=1 && o.padding<=32 &&
            o.maximum_spheres>0 && o.maximum_spheres<=static_cast<int>(ParticleGeometry::maximum_spheres),"Invalid GEMSS resolution/padding/sphere limit");
    require(std::isfinite(o.reconstruction_accuracy) && o.reconstruction_accuracy>0 && o.reconstruction_accuracy<=1 &&
            o.minimum_radius_voxels>=1 && o.minimum_radius_voxels<=o.voxel_resolution &&
            o.search_window>=1 && o.search_window<=o.voxel_resolution && o.persistence>=1 && o.persistence<=100 &&
            std::isfinite(o.minimum_center_distance_relative) && o.minimum_center_distance_relative>0 && o.minimum_center_distance_relative<=o.voxel_resolution &&
            std::isfinite(o.radius_offset_voxels) && o.radius_offset_voxels>=0 && o.radius_offset_voxels<=1,
            "Invalid GEMSS reconstruction parameters");
    require(std::isfinite(o.source_scale) && o.source_scale>0 && std::isfinite(o.density) && o.density>0 &&
            !o.source_units.empty() && o.maximum_voxels>0,"STL requires declared units, positive scale and density");
    const auto fingerprint=preflight(path);
    if (checksum) *checksum=fingerprint;
    auto mesh=GEMSS::load_mesh(path.string());
    for (int i=0;i<mesh.vertices.rows();++i) for (int j=0;j<3;++j) {
        const auto value=static_cast<long double>(mesh.vertices(i,j))*o.source_scale;
        require(std::isfinite(value) && std::abs(value)<=std::numeric_limits<float>::max(),"STL scale exceeds GEMSS float range");
        mesh.vertices(i,j)=static_cast<float>(value);
    }
    validate_mesh(mesh);
    const Eigen::Vector3f lo=mesh.vertices.colwise().minCoeff(),hi=mesh.vertices.colwise().maxCoeff();
    const Eigen::Vector3f extents=hi-lo;
    const float voxel=extents.minCoeff()/o.voxel_resolution;
    require(extents.allFinite() && extents.minCoeff()>0 && std::isfinite(voxel) && voxel>0,"Degenerate STL bounds or voxel size");
    require(mesh.vertices.cwiseAbs().maxCoeff()*std::numeric_limits<float>::epsilon()<voxel/16,
            "STL offset/scale loses voxel resolution in GEMSS float coordinates; recenter the input");
    long double count=1;
    for (int i=0;i<3;++i) count*=std::ceil(static_cast<long double>(extents[i])/voxel)+2*o.padding+1;
    require(count<=o.maximum_voxels && count< std::numeric_limits<int>::max(),"STL voxel grid exceeds configured resource limit");
    geometry_detail::arithmetic_range();
    return mesh;
}
}
void validate_stl(const std::filesystem::path& path,const GemssOptions& o) { checked_mesh(path,o); }
ParticleGeometry generate_stl(const std::filesystem::path& path,const GemssOptions& o) {
    geometry_detail::ValidationEnvironment environment;
    std::string checksum;
    auto mesh=checked_mesh(path,o,&checksum);
    const Eigen::Vector3f lo=mesh.vertices.colwise().minCoeff(),hi=mesh.vertices.colwise().maxCoeff();
    const Eigen::Vector3f extents=hi-lo;
    GEMSS::MultisphereConfig c;
    c.div=o.voxel_resolution; c.padding=o.padding; c.max_spheres=o.maximum_spheres;
    c.precision_target=o.reconstruction_accuracy; c.min_radius_vox=o.minimum_radius_voxels;
    c.search_window=o.search_window; c.persistence=o.persistence;
    c.min_center_distance_rel=o.minimum_center_distance_relative; c.radius_offset_vox=o.radius_offset_voxels;
    c.confine_mesh=o.confine_mesh; c.prune_isolated_spheres=o.prune_isolated_spheres;
    c.show_progress=false; c.is_2d=false; c.density=1; c.compute_physics=2; // Mandatory TARGET grid, never reconstruction union.
    // Use GEMSS's voxelization and reconstruction APIs. Keeping the target grid
    // also lets us reject empty masks and record the actual voxel spacing.
    geometry_detail::arithmetic_range();
    auto grid=GEMSS::mesh_to_binary_grid(mesh,c);
    require(grid.nz()>1 && std::any_of(grid.data.begin(),grid.data.end(),[](auto v){return v!=0;}),"No volumetric target voxels");
    auto pack=GEMSS::multisphere_from_voxels(grid,c);
    if (c.confine_mesh) {
        require(pack.centers.allFinite(),"Nonfinite GEMSS collision centers");
        // GEMSS's confinement wrapper uses libigl's default max<float> bound;
        // squaring that bound overflows even for a small valid cube. Use the
        // same libigl API with a finite geometry-sized search bound, and verify
        // distances before applying GEMSS's radius-minimum rule.
        const long double bound=2*std::hypot(static_cast<long double>(extents.x()),
            static_cast<long double>(extents.y()),static_cast<long double>(extents.z()));
        require(bound<std::sqrt(static_cast<long double>(std::numeric_limits<float>::max()))/4,
                "STL exceeds safe confinement distance range");
        Eigen::VectorXf distances; Eigen::VectorXi faces; Eigen::MatrixXf closest,normals;
        igl::signed_distance(pack.centers,mesh.vertices,mesh.triangles,igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL,
            0.0f,static_cast<float>(bound),distances,faces,closest,normals);
        require(distances.allFinite(),"Cannot establish finite mesh confinement distances");
        pack.radii=pack.radii.array().min(distances.array().abs());
    }
    require(pack.num_spheres()>0 && pack.num_spheres()<=static_cast<std::size_t>(o.maximum_spheres) &&
            std::isfinite(pack.precision) && pack.mass>0 && std::isfinite(pack.mass) && pack.center_of_mass.allFinite() && pack.inertia_tensor.allFinite(),
            "GEMSS failed to produce collision geometry and target physical properties");
    std::vector<Sphere> spheres;
    for (int i=0;i<pack.radii.size();++i)
        spheres.push_back({Vec3(pack.centers(i,0),pack.centers(i,1),pack.centers(i,2)),pack.radii[i]});
    // Pinned upstream samples at (i+0.5)*h, but physics sums i*h.
    // Restore the voxel-center coordinate convention shared by mesh and spheres.
    const long double half=static_cast<long double>(grid.voxel_size)/2;
    const Vec3 center(pack.center_of_mass.x()+half,pack.center_of_mass.y()+half,pack.center_of_mass.z()+half);
    PhysicalProperties p;
    p.volume=pack.mass; p.density=o.density;
    for (int i=0;i<3;++i) for (int j=0;j<3;++j) p.inertia_per_density[3*i+j]=pack.inertia_tensor(i,j);
    p.source=path.filename().string(); p.checksum=checksum; p.backend_version="GEMSS dc0f6de101b04545aa9679a13e692f60dfae2689";
    p.source_units=o.source_units; p.source_scale=o.source_scale; p.voxel_size=grid.voxel_size;
    std::ostringstream parameters; parameters.imbue(std::locale::classic()); parameters<<std::setprecision(9)
        <<"compute_physics=2;mesh_validation=closed_oriented_connected;self_intersections=unchecked;voxel_center_correction=+h/2;div="<<c.div<<";padding="<<c.padding
        <<";max_spheres="<<c.max_spheres<<";precision_target="<<c.precision_target<<";achieved="<<pack.precision
        <<";min_radius_vox="<<c.min_radius_vox<<";search_window="<<c.search_window<<";persistence="<<c.persistence
        <<";min_center_distance_rel="<<c.min_center_distance_rel<<";radius_offset_vox="<<c.radius_offset_vox
        <<";confine_mesh="<<c.confine_mesh<<";prune_isolated_spheres="<<c.prune_isolated_spheres;
    p.parameters=parameters.str();
    // Backend can set harmless underflow flags internally; validate its returned
    // finite values explicitly rather than exposing backend FP flags to the solver.
    require(!std::fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW),"GEMSS reported invalid arithmetic");
    std::feclearexcept(FE_ALL_EXCEPT);
    return ParticleGeometry(std::move(spheres),center,std::move(p),
        StlSource{std::filesystem::absolute(path).lexically_normal(),Vec3{},Quaternion{}},o.target_units);
}
}
