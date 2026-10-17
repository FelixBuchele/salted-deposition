// Adapted from the original SDBD source by Nikola Topic (2015).
// See Topic & Pöschel, JCP 308 (2016), DOI: 10.1016/j.jcp.2015.12.052.
#pragma once
#include <algorithm>
#include <cmath>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <gsl/gsl_rng.h>
#include "math/vector.hpp"
#include "math/periodic.hpp"
#include "math/matrix.hpp"
#include "math/rotation.hpp"
#include "invariant.hpp"
#include "log.hpp"
#include "sdbd/simulation.hpp"
namespace sdbd::detail {
class clustertype;
class surroundtype;
// One owner per run. Particles borrow this engine and cannot outlive the run.
class Engine {
public:
    // Create independent run state; copying is disabled because bodies borrow this owner.
    Engine() = default;
    // Copying is disabled because existing bodies borrow this owner.
    Engine(const Engine&) = delete;
    // Assignment cannot replace state borrowed by existing bodies.
    Engine& operator=(const Engine&) = delete;
    // Read the seven positional arguments and working-directory input, then deposit; return
    // zero or throw on invalid input/output.
    int run(int argc, char** argv);
    // Apply runtime settings and validate finite domains, supported grid spacing and RNG
    // name without drawing random numbers.
    void configure(const SimulationSettings&);
    // Validate grid dimensions, origin and equal integer spacings in legacy length units;
    // update periods or throw.
    void validate_parameters();
    // Validate particle/domain bounds and return imported geometry at dropping height;
    // update radius bound without RNG draws.
    std::optional<ParticleGeometry> prepare_particle(const std::optional<ParticleGeometry>&,int resolution,int count);
    // Initialize all bodies, then deposit sequentially using the seeded GSL stream; write
    // requested outputs and return zero, or throw.
    int execute(double requested_phi,int seedgsl,int resolution,int count,bool redrop,
                const std::string& identifier,const std::optional<ParticleGeometry>&,
                bool runtime_output=false,const std::function<void(int,int)>& progress={});

    // Draw an absolute horizontal COM position uniformly in area from a centered
    // disk of radius min(lx,ly)/4, consuming exactly two original GSL uniforms.
    vtype circular_drop_position();
    InjectionMode injection{InjectionMode::uniform};

    int N=500000;// particle count read from input, then overridden by the positional count

    int nx=200,ny=100,nz=1000;
    long double dx=100,dy=100,dz=100;
    long double xo=0,yo=0,zo=-2*dz;

    long double lx=nx*dx;
    long double ly=ny*dy;
    long double pxb=0;
    long double pxe=lx;
    long double pyb=0;
    long double pye=ly;

    std::string cmsfile="cms/cms";
    std::string positionsfile="positions/positions";
    std::string detailsfile="details/state";
    std::string particlefile; // Optional modern output: physical COM and retry counts.
    std::string vtkfile,stlposefile; // Empty disables optional final-state exports.
    std::string recordfile; // Consolidated run record; independent of numerical decisions.
    long double particle_scale{1},height_limit{};

    long double error=0;
    long double eps=60*1e-10;


    long double phi=0.01;// frustrated-motion angle; keep the original literal precision
    rotinftype rotinfphi;

    long double Rmax=34;//Maximal radius of a sphere in a system
    long double h0=60000000;//initial dropping height of particles
    long double epsz=60*1e-10;//particle is allowed to increase its height by maximum this value, to prevent issues when the structure is crystilline
    long double epstang=1e-6;//allows nearly grazing collisions to be classified as collisions.
    int particle{};

    long double l_pi=3.14159265359;


    // Compute downward travel dtl for two spheres (legacy length units); return false when
    // their horizontal disks do not intersect.
    bool pcontact0(vtype r1,long double R1,vtype r2,long double R2,long double& dtl);
    // Intersect a rotating sphere trajectory with a fixed sphere; return two
    // axis/cosine/sine solutions or false if no intersection exists.
    bool pcontact(vtype rp,const vtype& n,vtype rc,long double Rc,vtype rs,long double Rs,rotinftype& rotinf1,rotinftype& rotinf2);
    // Intersect circles of radii Rcd and Rsp separated by d; write planar points v1/v2 and
    // return whether strict intersections exist.
    bool circle_intersect(long double Rcd,long double Rsp,long double d,vtype& v1,vtype& v2);
    // Find two possible sphere-contact rotations about the horizontal steepest-descent axis
    // through rp; return false if absent.
    bool pcontact1(vtype rp,vtype cmo,vtype rcb,long double Rcb,vtype rsb,long double Rsb,rotinftype& rotinf1,rotinftype& rotinf2);
    // Write rotation bringing contacting constituent rc to support rs height, about the
    // COM-directed horizontal axis.
    void tohorizontal(const vtype& cm,const vtype& rc,const vtype& rs,rotinftype& rotinf);
    // Write rotation placing cm directly below support rp; require a defined horizontal
    // axis and positive sine.
    void tobelow1(const vtype& cm,const vtype& rp,rotinftype& rotinf);
    // Find rotations where constituent rc of radius Rc touches z=0 while rolling about rs;
    // write two solutions if present.
    bool toground1(const vtype&rs,vtype cm,vtype rc,const long double Rc,rotinftype& rotinf1,rotinftype& rotinf2);
    // Return dimensionless signed closing rate between rc and rs for one-support descent about rp; axis and distances must be nonzero.
    long double g(const vtype& rp, const vtype& cm, const vtype& rc, const vtype& rs);
    // Return the unit COM descent direction around support rp; require a nonzero horizontal offset.
    vtype grad(const vtype& rp,const vtype& cm);
    // Find sphere-contact rotations around distinct supports rp1/rp2, oriented toward
    // decreasing COM height; write two solutions if present.
    bool pcontact2(const vtype& rp1,const vtype& rp2,const vtype& cm,const vtype& rc,const long double Rcb,const vtype& rs,const long double Rsb,rotinftype& rotinf1,rotinftype& rotinf2);
    // Solve the grazing-contact condition (paper Eq. 25) for rotation on distinct supports;
    // write two rotations or return false.
    bool pdisconnect(vtype cm,vtype rp1,vtype rp2,vtype rd1,vtype rd2,rotinftype& rotinf1,rotinftype& rotinf2);
    // Find rotations placing rd1 at rp1 height about the rp1/rp2 axis; write solutions and
    // report numerical residual warnings.
    bool tohorizontal2(const vtype& rp1,const vtype& rp2,const vtype& cm,const vtype& rd1,rotinftype& rotinf1,rotinftype& rotinf2);
    // Write the hanging rotation placing cm below the two-support axis; require nonsingular
    // geometry.
    void tobelow2(const vtype& cm,const vtype& rp1,const vtype& rp2,rotinftype& rotinf);
    // Find rotations placing constituent rc of radius Rc on z=0 about two distinct
    // supports; write solutions if present.
    bool toground2(const vtype& rp1,vtype rp2,const vtype& cm, vtype rc,const long double Rc,rotinftype& rotinf1,rotinftype& rotinf2);
    // Return dimensionless signed closing rate of rc against rs for descent about two distinct supports.
    long double g2(const vtype& rp1,const vtype& rp2,const vtype& cm,const vtype& rc,const vtype& rs);
    // Return unit COM descent direction about two distinct supports; require a nonzero
    // rotational radius.
    vtype grad2(const vtype& rp1,const vtype& rp2,const vtype& cm);
    // Return whether initial/final COM offsets remain on the same side of the vertical
    // plane containing axis n.
    bool rightside(const vtype& cmr,const vtype& cmcr,const vtype n);
    // Intersect a radius-R circle normal to n with horizontal plane through rpl; write
    // sp/sm if real solutions exist (npl is unused).
    bool planecontact(vtype rp,vtype n,long double R,const vtype& rpl,vtype npl,vtype& sp,vtype& sm);
    // Write candidate grid index bounds for sphere r,R, including Rmax and horizontal
    // periodic images; equal grid spacings are required.
    void findindex(const vtype& r,long double R,int& ixl,int& ixr,int& iyl,int& iyr,int& izu);
    // Return rc distance in legacy length units from the horizontal descent axis through rp defined by cm.
    long double distance_from_axis(const vtype& rc,const vtype& rp,const vtype& cm);
    // Rotate r by phi radians about the horizontal descent axis through rp defined by cm;
    // return its new coordinates.
    vtype rotpoint(const vtype& r,const vtype& rp,const vtype& cm,long double phi);
    // Return rc distance in legacy length units from the axis through distinct points rp and rp2.
    long double distance_from_axis2(const vtype& rc,const vtype& rp,const vtype& rp2);
    // Rotate r by phi radians about the two-support descent axis oriented using cm; return
    // new coordinates.
    vtype rotpoint2(const vtype& r,const vtype& rp,const vtype& rp2,const vtype& cm,long double phi);
    // Own a newly allocated GSL generator of type, or the environment default when empty;
    // throw for invalid type/allocation.
    void initrand(const std::string& type = {});
    // Set the owned GSL generator seed without drawing; throw when the generator is uninitialized.
    void seed(unsigned long value);
    // Read the ordered legacy input fields including original output-label spelling;
    // validate or throw with field diagnostics.
    void load_parameters(const char* inputfile);
    // Return the sphere-volume-weighted COM of disjoint spheres r,R in input order; do not
    // use for overlapping-particle physics.
    vtype cmass(const std::vector<vtype>& r, const std::vector<long double>& R);
    // Check periodic interparticle overlaps using stored particle identity and legacy
    // tolerances; report/reject penetration, ignoring internal overlaps.
    int checkinters(surroundtype& su);
    // Check periodic interparticle overlaps using stored particle identity and legacy
    // tolerances; report/reject penetration, ignoring internal overlaps.
    bool checkinters(const clustertype& c,const surroundtype& s);
    // Append a settled body and its physical COM to the fixed environment; wrap constituent
    // coordinates into the periodic cell.
    void unload(clustertype& c,surroundtype& su);
    // Insert sphere v,R and its particle identity/active-boundary flag at index ii; reject
    // positions outside the finite grid.
    void unloadsphere(const vtype& v,int ii,long double R,surroundtype& su,std::size_t particle_id,bool collision_active);
    // Append unique fixed/moving coordinates and sphere-index pairs for active contacts within
    // eps_contact (legacy length units), preserving search order.
    void find_cluster_contacts(long double eps_contact,const clustertype& c,const surroundtype& su,std::vector<vtype>& rps,std::vector<vtype>& rds,std::vector<int>& p,std::vector<int>& d);
    // Write fixed sphere x,y,z,radius rows at legacy precision 15; preserve ordering and
    // throw on output failure.
    void output(surroundtype& su,std::string fname);
    // Append rigid-particle count and intersection status to name; throw on output failure.
    void writedetails(int size,int intersect,std::string name);
    // Append twice COM height z in legacy length units to name at precision 15; throw on
    // output failure.
    void write2cm(long double z,std::string name);
    // Draw one GSL uniform variate in [0,1); throw if the generator is uninitialized.
    long double ranf() { if (!rng_) throw std::logic_error("RNG not initialized"); return gsl_rng_uniform(rng_.get()); }
    // Return one signed uniform variate in [-1,1), consuming one owned GSL draw.
    long double ranfs() { return -1+2*ranf(); }
private:
    // Release the sole GSL generator when its owning run ends.
    struct RngDeleter {
        // Free rng through GSL; unique_ptr calls this only for a nonnull allocation.
        void operator()(gsl_rng* rng) const noexcept { gsl_rng_free(rng); }
    };
    std::unique_ptr<gsl_rng, RngDeleter> rng_;
};
}
