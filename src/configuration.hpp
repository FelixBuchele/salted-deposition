#pragma once
#include "sdbd/simulation.hpp"
#include "sdbd/gemss.hpp"
#include <map>

namespace sdbd::cli {
// One configuration value, with its file directory for resolving paths
// and its source line for error messages.
struct Value {
    std::string text;
    std::filesystem::path base;
    std::string location;
};
using Fields = std::map<std::string, Value>;
// STL input and generation settings, separate from deposition settings.
struct GenerationSettings {
    std::filesystem::path stl;
    GemssOptions options;
    long double scale{1};
};
// Settings after defaults, configuration and overrides have been combined.
// All lengths use units. Choose grid dimensions only after particle size is known.
struct Configuration {
    SimulationSettings simulation;
    std::string particle{"sphere"},units{"mm"},particle_checksum;
    std::optional<std::filesystem::path> particle_file;
    std::optional<GenerationSettings> generation;
    Vec3 box{720,720,6680};
    long double size{60},length{180},radius{30};
    std::size_t spheres{7};
    std::filesystem::path output;
    bool quiet{},legacy_output{},stl_poses{true};
};
// Read key=value fields and retain their source locations.
// Reject malformed lines and duplicate keys; paths are resolved later.
Fields read_configuration(const std::filesystem::path&);
// Apply defaults and interpret fields for a run or particle generation.
// Reject invalid settings or keys for the wrong operation. Output may be omitted
// for validation; this function does not generate geometry or draw random numbers.
Configuration resolve(const Fields&,bool generation,bool require_output=true);
// Derive equal integer grid spacings and an exact upper height from requested box;
// reject unsupported horizontal dimensions, resource limits or geometry mismatch.
void configure_domain(Configuration&,const std::optional<ParticleGeometry>&);
// Save effective inputs as a configuration that can be run again.
// Later metadata and final particle state are appended as comments.
void write_configuration(const std::filesystem::path&,const Configuration&);
}
