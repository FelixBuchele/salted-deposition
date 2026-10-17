#pragma once
#include "sdbd/particle.hpp"
namespace sdbd {
// Search grid: nx,ny,nz cells with spacings dx,dy,dz in simulation length units.
// The deposition plane is z=0; horizontal grid origins must be zero.
struct GridSettings {
    int nx{6},ny{6},nz{64};
    long double dx{120},dy{120},dz{120},xo{0},yo{0},zo{-1000};
};
// Initial horizontal placement: the original full-cell rule, or a physical COM
// sampled uniformly in a centered disk of radius min(Lx,Ly)/4.
enum class InjectionMode { uniform, circle };
// Settings for one run. Particles settle in sequence; accepted particles stay fixed.
struct SimulationSettings {
    GridSettings grid;
    double phi{0.01}; // Small rotation step in radians for frustrated contacts; original double precision.
    int seed{84328},resolution{1},shape{0},count{24};
    bool redrop{false};
    std::string rng; // Empty selects the GSL environment setting or default generator.
    std::function<void(int,int)> progress; // Called with completed and total particle counts.
    long double particle_scale{1}; // Built-in nominal dimensions relative to the original shape.
    long double height_limit{}; // Positive maximum accepted surface z; zero retains legacy bounds.
    InjectionMode injection{InjectionMode::uniform}; // Initial horizontal placement, including retries.

};
// Output paths; an empty path disables that file. The record path names an
// existing configuration to which final state and summary comments are appended.
struct SimulationFiles {
    std::filesystem::path positions,height,details,particles,vtk,stl_poses;
    std::filesystem::path record; // Run configuration receiving final COMs and retry counts.
};
// Resolve a GSL generator name from requested, the environment or the default.
// Reject unknown names without creating a generator or consuming random numbers.
std::string resolve_rng_name(std::string requested = {});
// Check settings and whether the optional particle fits the domain; throw on
// invalid input. This does not run deposition, draw random numbers or write files.
void validate_simulation(const SimulationSettings&,const std::optional<ParticleGeometry>& = std::nullopt);
// Run sequential deposition and write the requested files. Use the supplied
// geometry or the selected built-in shape; return zero or throw on failure.
int simulate(const SimulationSettings&,const SimulationFiles&,
             const std::optional<ParticleGeometry>& = std::nullopt);
// Read seven legacy arguments and working-directory input, then deposit and write original
// outputs; return zero or throw on invalid input or execution failure.
int simulate(int argc, char** argv);
}
