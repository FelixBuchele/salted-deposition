#pragma once
#include "sdbd/particle.hpp"
#include <span>

namespace sdbd {
// Write ASCII legacy PolyData points/vertices and radius, particle_id, sphere_id.
// Inputs are aligned final sphere arrays in solver units, grouped by increasing
// particle identity. Stream directly without an assembly copy; throw on bad data/I/O.
void write_deposition_vtk(const std::filesystem::path&,std::span<const Vec3> centers,
    std::span<const long double> radii,std::span<const std::size_t> particle_ids);
// Require retained GEMSS source frame, target-physics provenance and unchanged
// original STL bytes (FNV-1a64); throw rather than infer an absent source frame.
void validate_stl_source(const ParticleGeometry&);
// Return a deterministic FNV-1a64 fingerprint of file bytes; throw on read failure.
std::string file_checksum(const std::filesystem::path&);
// Write one point per sequential particle id: original STL origin P=C-R*c,
// normalized WXYZ orientation and uniform source scale. Centers use solver units.
// Orientations include the source-frame rotation. Caller retains the mesh reference
// in its run record; no mesh data or sidecar is written. Throw on invalid state/I/O.
void write_stl_poses(const std::filesystem::path&,std::span<const Vec3> centers,
    std::span<const Quaternion> orientations,const ParticleGeometry& source);
}
