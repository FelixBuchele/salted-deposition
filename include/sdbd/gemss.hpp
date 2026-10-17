#pragma once
#include "sdbd/particle.hpp"
namespace sdbd {
// Settings for the optional GEMSS backend. STL files do not declare units.
// source_scale converts input coordinates from source_units to target_units.
struct GemssOptions {
    int voxel_resolution{32}, padding{2}, maximum_spheres{64};
    float reconstruction_accuracy{0.95f};
    int minimum_radius_voxels{2}, search_window{2}, persistence{2};
    float minimum_center_distance_relative{1.4f}, radius_offset_voxels{0.5f};
    bool confine_mesh{false}, prune_isolated_spheres{false};
    long double source_scale{1}, density{1};
    std::string source_units;
    std::string target_units{"legacy"}; // Unit of generated coordinates after source_scale.
    std::size_t maximum_voxels{5'000'000};
};
// Generate collision spheres from binary STL using GEMSS and the declared units/scale.
// COM, volume and inertia come from the voxelized target mesh, not the sphere union.
// Return validated geometry or throw; never substitute a sphere-weighted COM.
ParticleGeometry generate_stl(const std::filesystem::path& source, const GemssOptions& options);
// Check parameters, binary STL topology, scaling and voxel resource bounds without
// voxelizing or generating spheres; throw on unsupported or invalid mesh/options.
void validate_stl(const std::filesystem::path& source,const GemssOptions& options);
}
