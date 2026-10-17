// Adapted from the original SDBD source by Nikola Topic (2015).
// See Topic & Pöschel, JCP 308 (2016), DOI: 10.1016/j.jcp.2015.12.052.
#include "engine.hpp"
#include <cstdlib>
#include <algorithm>
#include <numbers>
namespace {
// Find a named GSL generator type without allocating or drawing; throw when the registry has no match.
const gsl_rng_type* find_rng(const std::string& name) {
    for (auto types=gsl_rng_types_setup();*types;++types) if (name==(*types)->name) return *types;
    throw std::invalid_argument("Unknown GSL RNG type: "+name+"; use rng = mt19937 for the reference stream");
}
}
namespace sdbd {
std::string resolve_rng_name(std::string requested) {
    if (requested.empty()) {
        const auto* environment=std::getenv("GSL_RNG_TYPE");
        requested=environment ? environment : gsl_rng_mt19937->name;
    }
    return find_rng(requested)->name;
}
}
namespace sdbd::detail {
void Engine::initrand(const std::string& type) {
    const gsl_rng_type* selected;
    if (type.empty()) { gsl_rng_env_setup();selected=gsl_rng_default; }
    else selected=find_rng(type);
    rng_.reset(gsl_rng_alloc(selected));
    if (!rng_) throw std::bad_alloc();
}
vtype Engine::circular_drop_position() {
    const auto u_radius=ranf(),u_angle=ranf();
    const auto radius=std::min(lx,ly)/4*std::sqrt(u_radius);
    const auto angle=2*std::numbers::pi_v<long double>*u_angle;
    return vtype(lx/2+radius*std::cos(angle),ly/2+radius*std::sin(angle),0);
}
void Engine::seed(unsigned long value) {
    if (!rng_) throw std::logic_error("RNG not initialized");
    gsl_rng_set(rng_.get(), value);
}


}
