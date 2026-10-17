#pragma once
#include "sdbd/particle.hpp"
#include <span>

namespace sdbd {
// Write one VTK point per final sphere, with radius, particle_id and sphere_id.
// The three input arrays must match in length and be grouped by increasing particle id.
// Coordinates use simulation units. Write directly from the arrays without copying
// the assembly; throw on invalid data or file errors.
void write_deposition_vtk(const std::filesystem::path&,std::span<const Vec3> centers,
    std::span<const long double> radii,std::span<const std::size_t> particle_ids);
// Check that GEMSS target-based physical properties and the STL frame are present.
// Verify the original mesh against its stored checksum; throw if either check fails.
void validate_stl_source(const ParticleGeometry&);
// Return a deterministic FNV-1a64 fingerprint of file bytes; throw on read failure.
std::string file_checksum(const std::filesystem::path&);
// Write one VTK point per particle at the transformed STL origin, not at its COM.
// For final COM C, body-frame COM c and rotation R, the point is P=C-R*c, with c
// measured from the STL origin in simulation units. Include a unit (w,x,y,z)
// quaternion and uniform STL scale. Orientations already include the STL frame.
// The caller records the mesh path separately; no mesh is copied. Throw on invalid data/I/O.
void write_stl_poses(const std::filesystem::path&,std::span<const Vec3> centers,
    std::span<const Quaternion> orientations,const ParticleGeometry& source);
}
