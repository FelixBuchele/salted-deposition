// Built-in generator formulas and RNG call order retained from Nikola Topic's
// original 2015 SDBD implementation; see DOI: 10.1016/j.jcp.2015.12.052.
#include "sdbd/particle.hpp"
#include "geometry_validation.hpp"
#include <cmath>
#include <stdexcept>
namespace sdbd {
ParticleGeometry generate_stick(long double length, long double radius, std::size_t sphere_count) {
    geometry_detail::ValidationEnvironment environment;
    if (!std::isfinite(length) || !std::isfinite(radius) || radius<=0 ||
        length<=2*radius || sphere_count<2 || sphere_count>ParticleGeometry::maximum_spheres)
        throw std::invalid_argument("Stick requires finite length > 2*radius > 0 and 2..10000 spheres");
    const auto half=length/2-radius;
    const auto spacing=2*half/(sphere_count-1);
    if (!(spacing>0 && spacing<2*radius))
        throw std::invalid_argument("Stick spacing must be positive and less than the sphere diameter");
    std::vector<Sphere> spheres;
    spheres.reserve(sphere_count);
    // Symmetric index arithmetic gives paired centers with opposite coordinates.
    for (std::size_t i=0;i<sphere_count;++i)
        spheres.push_back({Vec3((2*static_cast<long double>(i)-(sphere_count-1))*half/(sphere_count-1),0,0),radius});
    geometry_detail::arithmetic_range();
    return ParticleGeometry(std::move(spheres),Vec3());
}
namespace {
using vtype = Vec3;
// Build the original shapes from sphere centers and radii.
// The caller supplies random numbers so geometry generation uses the simulation stream.
class LegacyGenerator {
public:
    std::vector<Vec3> r;
    std::vector<long double> R;
    const UniformRandom& random_;
    long double height_;
    long double scale_;
    // Return one validated uniform variate in [0,1), or throw for a bad callback result.
    long double ranf() {
        const auto value=random_();
        if (!std::isfinite(value) || value<0 || value>=1)
            throw std::invalid_argument("Generator random variate must be in [0,1)");
        return value;
    }
    // Return one signed uniform variate in [-1,1), consuming one callback draw.
    long double ranfs() { return -1+2*ranf(); }
    // Return a coordinate offset with components in [-amplitude,amplitude).
    // Draw z,y,x explicitly to keep the original random-number order.
    Vec3 jitter(long double amplitude) {
        const auto z=amplitude*ranfs();
        const auto y=amplitude*ranfs();
        const auto x=amplitude*ranfs();
        return Vec3(x,y,z);
    }
    // Create one sphere of radius 30*scale_ at height_, without random draws.
    void init_sphere();
    // Fill six perturbed axis arms at height_; nsph_per_side controls legacy discretization.
    void init_axes_lines_mshuffle(int nsph_per_side);
    // Fill perturbed cube surface lattice at height_, with nsph_per_side sites per edge.
    void init_cube_mshuffle(int nsph_per_side);
};
void LegacyGenerator::init_sphere()
{


	r.clear();R.clear();
	vtype v(0,0,height_);
	r.push_back(v);R.push_back(30.*scale_);


}

void LegacyGenerator::init_axes_lines_mshuffle(int nsph_per_side)
{
	r.clear();R.clear();
	vtype n;
	long double eps_sh=60*1e-4;
	long double Rp=(30.-eps_sh)*scale_;eps_sh*=scale_;
	long double le=2*Rp*(nsph_per_side-1);
	long double D=2*30.*scale_;

	n=1*ex;
	for(long double l=0;l<le+Rp/2.;l+=D) {vtype v=l*n;v.z()=v.z()+height_;const auto dv=jitter(eps_sh);r.push_back(v+dv);R.push_back(Rp-eps_sh*ranf());}

	n=-1*ex;
	for(long double l=D;l<le+Rp/2.;l+=D) {vtype v=l*n;v.z()=v.z()+height_;const auto dv=jitter(eps_sh);r.push_back(v+dv);R.push_back(Rp-eps_sh*ranf());}

	n=1*ey;
	for(long double l=D;l<le+Rp/2.;l+=D) {vtype v=l*n;v.z()=v.z()+height_;const auto dv=jitter(eps_sh);r.push_back(v+dv);R.push_back(Rp-eps_sh*ranf());}

	n=-1*ey;
	for(long double l=D;l<le+Rp/2.;l+=D) {vtype v=l*n;v.z()=v.z()+height_;const auto dv=jitter(eps_sh);r.push_back(v+dv);R.push_back(Rp-eps_sh*ranf());}

	n=1*ez;
	for(long double l=D;l<le+Rp/2.;l+=D) {vtype v=l*n;v.z()=v.z()+height_;const auto dv=jitter(eps_sh);r.push_back(v+dv);R.push_back(Rp-eps_sh*ranf());}

	n=-1*ez;
	for(long double l=D;l<le+Rp/2.;l+=D) {vtype v=l*n;v.z()=v.z()+height_;const auto dv=jitter(eps_sh);r.push_back(v+dv);R.push_back(Rp-eps_sh*ranf());}

}

void LegacyGenerator::init_cube_mshuffle(int nsph_per_side)
{
	r.clear();R.clear();
	vtype n;
	long double eps_sh=60*1e-4;
	long double Rp=(30.-eps_sh)*scale_;eps_sh*=scale_;
	long double D=2*30.*scale_;

	// Bottom face: include every edge and corner.
	for(int i=0;i<nsph_per_side;i++)
	for(int j=0;j<nsph_per_side;j++)
	{vtype v=D*(i*ex+j*ey);v.z()=v.z()+height_;const auto dv=jitter(eps_sh);r.push_back(v+dv);R.push_back(Rp-eps_sh*ranf());}

	// Top face: same lattice, shifted by the cube edge length.
	for(int i=0;i<nsph_per_side;i++)
	for(int j=0;j<nsph_per_side;j++)
	{vtype v=D*(i*ex+j*ey);v.z()=v.z()+height_;v.z()=v.z()+(nsph_per_side-1)*D;const auto dv=jitter(eps_sh);r.push_back(v+dv);R.push_back(Rp-eps_sh*ranf());}

	// Front face: skip edges already included on the bottom and top.
	for(int i=0;i<nsph_per_side;i++)
	for(int j=1;j<nsph_per_side-1;j++)
	{vtype v=D*(i*ex+j*ez);v.z()=v.z()+height_;const auto dv=jitter(eps_sh);r.push_back(v+dv);R.push_back(Rp-eps_sh*ranf());}

	// Back face: the same edge exclusion as the front.
	for(int i=0;i<nsph_per_side;i++)
	for(int j=1;j<nsph_per_side-1;j++)
	{vtype v=D*(i*ex+j*ez);v.z()=v.z()+height_;v.y()=v.y()+(nsph_per_side-1)*D;const auto dv=jitter(eps_sh);r.push_back(v+dv);R.push_back(Rp-eps_sh*ranf());}

	// Left face: skip edges already included by the other four faces.
	for(int i=1;i<nsph_per_side-1;i++)
	for(int j=1;j<nsph_per_side-1;j++)
	{vtype v=D*(i*ey+j*ez);v.z()=v.z()+height_;const auto dv=jitter(eps_sh);r.push_back(v+dv);R.push_back(Rp-eps_sh*ranf());}

	// Right face: likewise include only the remaining face interior.
	for(int i=1;i<nsph_per_side-1;i++)
	for(int j=1;j<nsph_per_side-1;j++)
	{vtype v=D*(i*ey+j*ez);v.z()=v.z()+height_;v.x()=v.x()+(nsph_per_side-1)*D;const auto dv=jitter(eps_sh);r.push_back(v+dv);R.push_back(Rp-eps_sh*ranf());}

}

}
ParticleGeometry generate_builtin(BuiltinShape shape, int resolution,
                                  const UniformRandom& random, long double initial_height,long double length_scale) {
    geometry_detail::ValidationEnvironment environment;
    if (resolution<=0 || !std::isfinite(initial_height) || !std::isfinite(length_scale) || length_scale<=0 ||
        (shape==BuiltinShape::cube && resolution<2))
        throw std::invalid_argument("Invalid built-in resolution or initial height");
    const long double n=resolution;
    const long double count=shape==BuiltinShape::cube ? 6*n*n-12*n+8 :
                            shape==BuiltinShape::spiky ? 6*n-5 : 1;
    if (count>ParticleGeometry::maximum_spheres)
        throw std::invalid_argument("Built-in sphere count exceeds supported limit");
    if (shape!=BuiltinShape::sphere && !random)
        throw std::invalid_argument("Built-in generator requires a random source");
    LegacyGenerator generator{{},{},random,initial_height,length_scale};
    switch (shape) {
        case BuiltinShape::sphere: generator.init_sphere(); break;
        case BuiltinShape::spiky: generator.init_axes_lines_mshuffle(resolution); break;
        case BuiltinShape::cube: generator.init_cube_mshuffle(resolution); break;
        default: throw std::invalid_argument("Unknown built-in shape");
    }
    std::vector<Sphere> spheres;
    spheres.reserve(generator.r.size());
    for (std::size_t i=0;i<generator.r.size();++i)
        spheres.push_back({generator.r[i],generator.R[i]});
    geometry_detail::arithmetic_range();
    return ParticleGeometry(std::move(spheres));
}
}
