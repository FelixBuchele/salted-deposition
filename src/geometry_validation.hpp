#pragma once
#include <cfenv>
#include <stdexcept>
namespace sdbd::geometry_detail {
// Geometry import/validation must report range errors even when the caller has
// enabled solver FP traps. Restore both control settings and flags on exit.
class ValidationEnvironment {
    std::fenv_t saved_{};
public:
    // Mask floating-point traps and clear flags during validation; throw if the environment
    // cannot be saved.
    ValidationEnvironment() {
        if (std::feholdexcept(&saved_)!=0) throw std::runtime_error("Cannot validate floating-point geometry");
    }
    // Restore the caller's saved floating-point controls and flags.
    ~ValidationEnvironment() { std::fesetenv(&saved_); }
    // Environment guards cannot be copied; each instance restores its own saved state.
    ValidationEnvironment(const ValidationEnvironment&)=delete;
    // Assignment is disabled to prevent replacing a live environment guard.
    ValidationEnvironment& operator=(const ValidationEnvironment&)=delete;
};
// Reject invalid/divide/overflow/underflow flags accumulated during geometry arithmetic.
inline void arithmetic_range() {
    if (std::fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW))
        throw std::invalid_argument("Geometry arithmetic exceeds numerical range");
}
}
