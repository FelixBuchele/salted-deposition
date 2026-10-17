// Adapted from the original SDBD source by Nikola Topic (2015).
// See Topic & Pöschel, JCP 308 (2016), DOI: 10.1016/j.jcp.2015.12.052.
#pragma once
#include <cmath>
#include <iostream>
#include "vector.hpp"

namespace sdbd::detail {
// A 3x3 matrix used to convert local circle coordinates to simulation coordinates.
class matrix
{
	public:
	// Construct a 3x3 coordinate transform whose columns are v1,v2,v3; store long-double components.
	matrix(const vtype& v1,const vtype& v2,const vtype& v3);

	public:
	long double m11,m12,m13;
	long double m21,m22,m23;
	long double m31,m32,m33;
};

// Return matrix m applied to vector v using the original multiplication and summation order.
vtype ma(const matrix& m, const vtype& v);
}
