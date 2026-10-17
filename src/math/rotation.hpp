// Adapted from the original SDBD source by Nikola Topic (2015).
// See Topic & Pöschel, JCP 308 (2016), DOI: 10.1016/j.jcp.2015.12.052.
#pragma once
#include <cmath>
#include <iostream>
#include "vector.hpp"

namespace sdbd::detail {
// Axis-angle rotation stored as unit axis n and angle cosine cc/sine ss.
class rotinftype
{
	public:
	long double cc{},ss{};
	vtype n;
};

// Compose rotations about the same axis using sine/cosine addition; reject differing axes.
rotinftype rotn(const rotinftype& rotinf1,const rotinftype& rotinf2);
}
