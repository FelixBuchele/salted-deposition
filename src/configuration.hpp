#pragma once
#include "sdbd/simulation.hpp"
#include "sdbd/gemss.hpp"
#include <map>

namespace sdbd::cli {
// Parsed value with file-relative paths and line context for diagnostics.
struct Value {
    std::string text;
    std::filesystem::path base;
    std::string location;
};
using Fields = std::map<std::string, Value>;
// Geometry preparation only; source_scale converts source units into output units.
struct GenerationSettings {
    std::filesystem::path stl;
    GemssOptions options;
    long double scale{1};
};
// One resolved run or generation configuration; box dimensions use units throughout.
// Grid derivation is deferred until the prepared geometry's bounds are known.
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
// Read key=value configuration without resolving paths; reject malformed/duplicate fields.
Fields read_configuration(const std::filesystem::path&);
// Resolve run/generation keys and units, rejecting invalid or misplaced fields.
// Validation may omit output; no mesh generation or random draws occur here.
Configuration resolve(const Fields&,bool generation,bool require_output=true);
// Derive equal integer grid spacings and an exact upper height from requested box;
// reject unsupported horizontal dimensions, resource limits or geometry mismatch.
void configure_domain(Configuration&,const std::optional<ParticleGeometry>&);
// Save replayable inputs; metadata and deposited state are appended as ignored comments.
void write_configuration(const std::filesystem::path&,const Configuration&);
}
