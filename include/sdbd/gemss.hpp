#pragma once
#include "sdbd/particle.hpp"
namespace sdbd {
// Available only from the optional sdbd_gemss target. STL has no intrinsic units;
// source_scale converts declared source_units to the solver's legacy length unit.
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
// Generate collision spheres from binary STL via pinned GEMSS; options declare units/scaling.
// Return target-voxel COM/properties or throw, never sphere-weight fallback.
ParticleGeometry generate_stl(const std::filesystem::path& source, const GemssOptions& options);
// Check parameters, binary STL topology, scaling and voxel resource bounds without
// voxelizing or generating spheres; throw on unsupported or invalid mesh/options.
void validate_stl(const std::filesystem::path& source,const GemssOptions& options);
}
