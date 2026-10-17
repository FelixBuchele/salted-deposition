// Adapted from the original SDBD source by Nikola Topic (2015).
// See Topic & Pöschel, JCP 308 (2016), DOI: 10.1016/j.jcp.2015.12.052.
#pragma once
#include <cmath>
#include <iostream>
#include "vector.hpp"

namespace sdbd::detail {
//implements modular operations needed for periodic boundary conditions

// Return i modulo positive n in [0,n), including negative indices for periodic grid lookup.
int imod(int i,int n);
// Return x modulo positive period l in [0,l), with long-double remainder in the original expression order.
long double rmod(long double x,long double l);
// Return nearest periodic displacement in [-l/2,l/2) for positive period l.
long double moddist(long double x, long double l);
// Return horizontal periodic wrapping, or the nearest image relative to rb; preserve z and
// require positive periods lx,ly.
vtype vmod(const vtype& r,long double lx,long double ly);
// Return horizontal periodic wrapping, or the nearest image relative to rb; preserve z and
// require positive periods lx,ly.
vtype vmod(const vtype& ri,const vtype& rb,long double lx,long double ly);
}
