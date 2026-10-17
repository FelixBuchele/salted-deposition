#pragma once
#include <stdexcept>
// Remains active in Release builds; preserves the legacy assertion predicate.
#define SDBD_REQUIRE(condition) do { if (!(condition)) throw std::runtime_error("Invariant failed: " #condition); } while (false)
