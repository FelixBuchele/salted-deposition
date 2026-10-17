#pragma once
#include "sdbd/particle.hpp"
namespace sdbd {
// Finite search-grid dimensions/spacings in legacy length units; deposition plane is z=0 and
// horizontal origins must be zero.
struct GridSettings {
    int nx{6},ny{6},nz{64};
    long double dx{120},dy{120},dz{120},xo{0},yo{0},zo{-1000};
};
// Initial horizontal COM sampling: full periodic cell or centered disk of radius min(Lx,Ly)/4.
enum class InjectionMode { uniform, circle };
// Runtime deposition parameters; preserve GSL stream selection, double-parsed phi and
// immutable-deposition model.
struct SimulationSettings {
    GridSettings grid;
    double phi{0.01}; // Parse as double, then widen, matching the original input.
    int seed{84328},resolution{1},shape{0},count{24};
    bool redrop{false};
    std::string rng; // Empty retains environment-based selection for API compatibility.
    std::function<void(int,int)> progress; // Completed particles, total; no RNG draws.
    long double particle_scale{1}; // Built-in nominal dimensions relative to the original shape.
    long double height_limit{}; // Positive maximum accepted surface z; zero retains legacy bounds.
    InjectionMode injection{InjectionMode::uniform}; // Initial horizontal placement, including retries.

};
// Requested outputs; empty paths disable individual files. record appends final
// state to an existing configuration, without affecting numerical routines.
struct SimulationFiles {
    std::filesystem::path positions,height,details,particles,vtk,stl_poses;
    std::filesystem::path record; // Optional existing run.conf to append final COM/retry/summary comments.
};
// Return a validated GSL RNG name from requested or environment/default, without allocating or
// drawing; throw for unknown types.
std::string resolve_rng_name(std::string requested = {});
// Validate settings and optional geometry/domain compatibility without RNG draws, deposition or
// output writes; throw on invalid input.
void validate_simulation(const SimulationSettings&,const std::optional<ParticleGeometry>& = std::nullopt);
// Deposit sequentially with supplied settings/geometry and files, or seven legacy arguments
// plus working-directory input; return zero or throw on failure.
int simulate(const SimulationSettings&,const SimulationFiles&,
             const std::optional<ParticleGeometry>& = std::nullopt);
// Read seven legacy arguments and working-directory input, then deposit and write original
// outputs; return zero or throw on invalid input or execution failure.
int simulate(int argc, char** argv);
}
