// Adapted from the original SDBD source by Nikola Topic (2015).
// See Topic & Pöschel, JCP 308 (2016), DOI: 10.1016/j.jcp.2015.12.052.
#pragma once
#include "math/vector.hpp"
#include "sdbd/orientation.hpp"
#include <cstddef>
#include <vector>
namespace sdbd::detail {
class Engine;
// Fixed deposited spheres, particle identities, physical COMs and spatial search grid; all
// sphere-indexed arrays remain aligned.
class surroundtype
{
	public:
	std::vector<vtype> r;
	std::vector<long double> R;
    std::vector<std::size_t> particle_ids;
    std::vector<bool> collision_active;
    std::vector<vtype> particle_centers;
    std::vector<Quaternion> particle_orientations; // Only accepted poses when requested.
	std::vector<std::vector<std::vector<std::vector<int> > > > box;
	std::vector<std::vector<int> > izmax;
	// Allocate the fixed-environment grid and empty-column height indices from engine dimensions;
	// propagate allocation failures.
	void init_grid(const Engine& engine);
};
}
