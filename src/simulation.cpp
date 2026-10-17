// Adapted from the original SDBD source by Nikola Topic (2015).
// See Topic & Pöschel, JCP 308 (2016), DOI: 10.1016/j.jcp.2015.12.052.
#include "sdbd/simulation.hpp"
#include "sdbd/export.hpp"
#include "engine.hpp"
#include "cluster.hpp"
#include "surround.hpp"
#include "geometry_validation.hpp"
#include <fenv.h>
#include <utility>
#include <cerrno>
#include <climits>
#include <cstdlib>
#include <limits>
#include <string_view>
#include <optional>
#include <iomanip>
#include <locale>
namespace sdbd::detail {
namespace {
// Parse the complete positional integer text within int range; throw with parameter name on failure.
int integer(const char* text, const char* name) {
    errno = 0;
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (errno == ERANGE || end == text || *end != '\0' || value < INT_MIN || value > INT_MAX)
        throw std::invalid_argument(std::string("Invalid integer: ") + name);
    return static_cast<int>(value);
}
// Parse a positive finite angle in radians as double, preserving legacy input precision; throw otherwise.
double angle(const char* text) {
    errno = 0;
    char* end = nullptr;
    // Deliberately parse as double, then widen: the legacy program used atof.
    const double value = std::strtod(text, &end);
    if (errno == ERANGE || end == text || *end != '\0' || !std::isfinite(value) || value <= 0)
        throw std::invalid_argument("Invalid positive finite phi");
    return value;
}
// Restore caller floating-point state on success and exception unwinding.
class FloatingEnvironment {
    fenv_t saved_{};
public:
    // Save caller state and enable original invalid/divide/overflow/underflow traps; throw
    // if unavailable.
    FloatingEnvironment() {
        if (fegetenv(&saved_) != 0) throw std::runtime_error("Cannot save floating-point environment");
        if (feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW) == -1)
            throw std::runtime_error("Cannot enable floating-point traps");
    }
    // Restore caller floating-point state on normal exit or exception unwinding.
    ~FloatingEnvironment() { fesetenv(&saved_); }
    // Environment guards cannot be copied; each instance restores its own saved state.
    FloatingEnvironment(const FloatingEnvironment&) = delete;
    // Assignment is disabled to prevent replacing a live environment guard.
    FloatingEnvironment& operator=(const FloatingEnvironment&) = delete;
};
}
int Engine::run(int argc,char** argv)
{
    if (argc != 8) throw std::invalid_argument("Usage: deposit phi seed resolution shape count redrop identifier");
    const double requested_phi = angle(argv[1]);
    const int seedgsl = integer(argv[2], "seed");
    const int nspheres = integer(argv[3], "resolution");
    const std::string_view shape_argument=argv[4];
    const bool manual=shape_argument.starts_with("file:");
    particle = manual ? 3 : integer(argv[4], "shape");
    const int count = integer(argv[5], "count");
    const int redrop = integer(argv[6], "redrop");
    if (nspheres <= 0 || count <= 0 || particle < 0 || (!manual && particle > 2) ||
        (particle == 2 && nspheres < 2) || (redrop != 0 && redrop != 1))
        throw std::invalid_argument("Invalid resolution, shape, count or redrop setting");
    const bool doredrop = redrop != 0;
    const std::string file_identifier = argv[7];
    initrand();
    load_parameters("input");
    std::optional<ParticleGeometry> geometry;
    if (manual) {
        if (shape_argument.size()==5) throw std::invalid_argument("Missing particle file path");
        geometry=load_particle(std::filesystem::path(shape_argument.substr(5)));
    }
    return execute(requested_phi,seedgsl,nspheres,count,doredrop,file_identifier,geometry);
}

void Engine::configure(const SimulationSettings& settings) {
    geometry_detail::ValidationEnvironment environment;
    if (settings.injection!=InjectionMode::uniform && settings.injection!=InjectionMode::circle)
        throw std::invalid_argument("Unknown injection mode");
    injection=settings.injection;
    N=settings.count;phi=settings.phi;
    nx=settings.grid.nx;ny=settings.grid.ny;nz=settings.grid.nz;
    dx=settings.grid.dx;dy=settings.grid.dy;dz=settings.grid.dz;
    xo=settings.grid.xo;yo=settings.grid.yo;zo=settings.grid.zo;
    if (!std::isfinite(settings.phi) || settings.phi<=0 || settings.resolution<=0 ||
        settings.shape<0 || settings.shape>2 || (settings.shape==2 && settings.resolution<2))
        throw std::invalid_argument("Invalid phi, shape or resolution");
    validate_parameters();
    if (std::numeric_limits<long double>::epsilon()*std::max({lx,ly,h0})>eps)
        throw std::invalid_argument("Domain coordinates cannot resolve the solver tolerance");
    if (!settings.rng.empty()) resolve_rng_name(settings.rng);
    particle=settings.shape;
    if (!std::isfinite(settings.particle_scale) || settings.particle_scale<=0 ||
        !std::isfinite(settings.height_limit) || settings.height_limit<0)
        throw std::invalid_argument("Invalid particle scale or maximum surface height");
    particle_scale=settings.particle_scale;height_limit=settings.height_limit;
    Rmax=std::max(Rmax,34*particle_scale);
    geometry_detail::arithmetic_range();
}

std::optional<ParticleGeometry> Engine::prepare_particle(const std::optional<ParticleGeometry>& input_geometry,int nspheres,int count) {
    geometry_detail::ValidationEnvironment environment;
    std::optional<ParticleGeometry> manual_geometry;
    if (input_geometry) {
        const auto& geometry=*input_geometry;
        Rmax=std::max(Rmax,1.1L*geometry.maximum_radius());
        const auto extent=geometry.bounding_radius();
        if (!std::isfinite(Rmax) || 2*extent>=std::min(lx,ly) ||
            2*extent>=zo+nz*dz || 2*extent>=h0 || -zo<=2*extent+2*Rmax ||
            std::numeric_limits<long double>::epsilon()*std::max({lx,ly,h0})>eps)
            throw std::invalid_argument("Manual particle exceeds domain bounds or below-ground search margin");
        std::vector<Sphere> placed;
        for (const auto& sphere : geometry.spheres())
            placed.push_back({sphere.center-geometry.center_of_mass()+vtype(0,0,h0),sphere.radius});
        auto source=geometry.stl_source();
        if (source) source->origin+=vtype(0,0,h0)-geometry.center_of_mass();
        manual_geometry.emplace(std::move(placed),vtype(0,0,h0),geometry.physical_properties(),std::move(source));
        // Match the original centered-sphere handling for an exactly spherical
        // union, including redundant contained spheres, without deleting geometry.
        if (geometry.is_centered_sphere()) particle=0;
    }
    if (!input_geometry) {
        const long double extent=particle==0 ? 0 : 60.*(nspheres-1)*particle_scale;
        if (!std::isfinite(extent) || 2*extent+60*particle_scale>=std::min(lx,ly))
            throw std::invalid_argument("Particle resolution exceeds horizontal domain bounds");
    }
    long double spheres=1;
    if (input_geometry) spheres=input_geometry->spheres().size();
    else if (particle==1) spheres=6.L*(nspheres-1)+1;
    else if (particle==2) spheres=6.L*nspheres*nspheres-12.L*nspheres+8;
    if (spheres>ParticleGeometry::maximum_spheres || spheres*count>INT_MAX)
        throw std::invalid_argument("Particle sphere count exceeds supported resource/index range");
    geometry_detail::arithmetic_range();
    return manual_geometry;
}

int Engine::execute(double requested_phi,int seedgsl,int nspheres,int count,bool doredrop,
                    const std::string& file_identifier,const std::optional<ParticleGeometry>& geometry,
                    bool runtime_output,const std::function<void(int,int)>& progress) {
    const bool manual=geometry.has_value();
    if (manual) particle=3;
    const auto manual_geometry=prepare_particle(geometry,nspheres,count);
    FloatingEnvironment floating_environment;
    phi = requested_phi;
    rotinfphi.cc=std::cos(phi);
    rotinfphi.ss=std::sin(phi);
    if (!runtime_output) ds("phi",phi);
    seed(seedgsl);
    N = count;

	if (!runtime_output) {
	if(manual) ds("Manual multisphere particles");
	else if(particle==0) ds("Spheres");
	else if(particle==1) ds("Spiky particles");
	else if(particle==2) ds("Cubes");
    }

	std::vector<clustertype> co;
	co.reserve(N);

	for(int i=0;i<N;i++) //filling up the vector containing clusters, represented with object cp. The particles from this vector will be deposited
	{
        const auto generated=manual ? std::optional<ParticleGeometry>{} : std::optional<ParticleGeometry>{generate_builtin(
            static_cast<BuiltinShape>(particle),nspheres,[this]{ return ranf(); },h0,particle_scale)};
        const auto& definition=manual ? *manual_geometry : *generated;
        if (definition.spheres().size()>static_cast<std::size_t>(INT_MAX)/static_cast<std::size_t>(N))
            throw std::invalid_argument("Total sphere count exceeds solver index range");
        clustertype cp(*this,definition,static_cast<std::size_t>(i));

        long double ddx,ddy;
        if (injection==InjectionMode::circle) {
            const auto point=circular_drop_position();
            // Center the physical COM; imported STL origins need not be centered.
            ddx=point.x()-cp.cm.x();ddy=point.y()-cp.cm.y();
        } else {
            // Retain the original expressions and draw order for uniform injection.
            ddx=lx/2+lx/2*ranfs();ddy=ly/2+ly/2*ranfs();
        }
        cp.shift_horizontal(ddx,ddy,!manual);
		cp.randrot();

		co.push_back(std::move(cp));	
	}

	surroundtype su;
	su.init_grid(*this);//initializations of the object containing fixed particles
    std::vector<std::size_t> redrops(particlefile.empty() && recordfile.empty() ? 0 : N,0);

	for(int i=0;i<co.size();i++)//loop in which particles are sedimented 
	{
		if(!runtime_output && i%100==0) ds("i",i);
		clustertype c=std::move(co[i]);
		if(!doredrop)//if no redropping leave particle in (rare) unstable positions 
		{
			c.settle(su);//sediments the particle (cluster)		
		}
		else//if there is a problem redropping of a cluster is allowed
		{
			bool collapsed;
			clustertype cplr=c;
			do//loop until the particle reaches a proper stable position (needed very rarely)
			{
				collapsed=c.settle(su);//sediments the particle (cluster)		
				if(collapsed)//choose new random position and orientation
				{
                    if (!redrops.empty()) ++redrops[i];
					c=cplr;
                    long double ddx,ddy;
                    if (injection==InjectionMode::circle) {
                        const auto point=circular_drop_position();
                        ddx=point.x()-c.cm.x();ddy=point.y()-c.cm.y();
                    } else {
                        ddx=-c.cm.x()+lx/2.+lx/2.*ranfs();
                        ddy=-c.cm.y()+ly/2.+ly/2.*ranfs();
                    }

                    c.shift_horizontal(ddx,ddy,!manual);
					c.randrot();	
				}
			}
			while(collapsed);
		}
        if (height_limit>0) {
            for (std::size_t j=0;j<c.r.size();++j)
                if (c.r[j].z()+c.R[j]>height_limit)
                    throw std::out_of_range("Deposited particle surface exceeds requested box height");
        }
		unload(c,su);//the cluster has been sedimented, add it to the fixed particles
        if (progress) progress(i+1,N);
	}

	//output information about the simulation

	std::string positionsfile_str=positionsfile+file_identifier;
	if (!positionsfile.empty()) output(su,positionsfile_str.c_str());//positions

	std::string cmsfile_str=cmsfile+file_identifier;
    long double packing_center_z;
    if (manual) {
        // Identical copies of one file geometry have equal physical mass. Their
        // explicit COMs define the packing COM even if constituent spheres overlap.
        long double mean=0;
        for (const auto& center : su.particle_centers) mean+=center.z();
        packing_center_z=mean/su.particle_centers.size();
    } else packing_center_z=cmass(su.r,su.R).z();
    if (!cmsfile.empty()) write2cm(packing_center_z,cmsfile_str.c_str());//filling height

	int intersect=checkinters(su);//check for intersections
	std::string detailsfile_str=detailsfile+file_identifier;
	if (!detailsfile.empty()) writedetails(co.size(),intersect,detailsfile_str.c_str());//other information
    if (!particlefile.empty()) {
        std::ofstream file(particlefile);
        file.imbue(std::locale::classic());
        file<<std::setprecision(std::numeric_limits<long double>::max_digits10);
        file<<"# particle_id com_x com_y com_z redrops; coordinates in legacy length units\n";
        for (std::size_t i=0;i<su.particle_centers.size();++i)
            file<<i<<' '<<su.particle_centers[i]<<' '<<redrops[i]<<'\n';
        file.close();
        if (!file) throw std::runtime_error("Cannot write physical particle states: "+particlefile);
    }

    if (!vtkfile.empty()) write_deposition_vtk(vtkfile,su.r,su.R,su.particle_ids);
    if (!stlposefile.empty()) write_stl_poses(stlposefile,su.particle_centers,su.particle_orientations,*geometry);
    if (!recordfile.empty()) {
        std::ofstream file(recordfile,std::ios::app);file.imbue(std::locale::classic());
        file<<std::setprecision(std::numeric_limits<long double>::max_digits10);
        long double top=0;
        for (std::size_t i=0;i<su.r.size();++i) top=std::max(top,su.r[i].z()+su.R[i]);
        file<<"# surface_height "<<top<<"\n# twice_packing_com_z "<<2*packing_center_z
            <<"\n# intersections "<<intersect<<"\n# particle_id com_x com_y com_z redrops\n";
        for (std::size_t i=0;i<su.particle_centers.size();++i)
            file<<"# particle "<<i<<' '<<su.particle_centers[i]<<' '<<redrops[i]<<'\n';
        file.close();if (!file) throw std::runtime_error("Cannot append deposited state to run record");
    }
	return 0;
}
}
namespace sdbd {
void validate_simulation(const SimulationSettings& settings,const std::optional<ParticleGeometry>& geometry) {
    detail::Engine engine;engine.configure(settings);
    if (geometry) engine.particle=3;
    engine.prepare_particle(geometry,settings.resolution,settings.count);
}
int simulate(const SimulationSettings& settings,const SimulationFiles& files,const std::optional<ParticleGeometry>& geometry) {
    detail::Engine engine;engine.configure(settings);
    engine.positionsfile=files.positions.string();engine.cmsfile=files.height.string();engine.detailsfile=files.details.string();
    engine.particlefile=files.particles.string();
    engine.vtkfile=files.vtk.string();engine.stlposefile=files.stl_poses.string();
    engine.recordfile=files.record.string();
    if (!files.stl_poses.empty()) {
        if (!geometry) throw std::invalid_argument("STL pose export is unavailable for built-in geometry");
        validate_stl_source(*geometry);
    }
    engine.initrand(settings.rng);
    return engine.execute(settings.phi,settings.seed,settings.resolution,settings.count,settings.redrop,"",geometry,true,settings.progress);
}
int simulate(int argc, char** argv) { detail::Engine engine; return engine.run(argc, argv); }
}
