#pragma once
#include "sdbd/vector.hpp"
#include <array>

namespace sdbd {
// Unit quaternion in scalar-first order (w,x,y,z), using right-handed rotations.
// It rotates vectors in fixed axes. The solver tracks it for export, but continues
// to calculate motion with its original rotation formulas.
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
    // Convert an axis and angle cosine/sine into a unit quaternion; normalize the axis.
    static Quaternion from_axis(Vec3 axis,long double cosine,long double sine);
};
}
