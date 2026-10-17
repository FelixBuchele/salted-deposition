#pragma once
#include "sdbd/vector.hpp"
#include <array>

namespace sdbd {
// Unit right-handed active rotation, scalar-first WXYZ; independent export bookkeeping.
struct Quaternion {
    long double w{1},x{},y{},z{};
    // Return a unit quaternion or throw for nonfinite/zero norm.
    Quaternion normalized() const;
    // Rotate a length vector by this unit quaternion; does not translate it.
    Vec3 apply(Vec3 value) const;
    // Return the inverse active rotation of a unit quaternion.
    Quaternion inverse() const { return {w,-x,-y,-z}; }
    // Return the unit composition a*b, applying b before a.
    friend Quaternion operator*(const Quaternion& a,const Quaternion& b);
    // Convert orthonormal rotation columns into a normalized WXYZ quaternion.
    static Quaternion from_basis(const std::array<Vec3,3>& columns);
    // Convert a unit-axis rotation's sine/cosine into a normalized quaternion.
    static Quaternion from_axis(Vec3 axis,long double cosine,long double sine);
};
}
