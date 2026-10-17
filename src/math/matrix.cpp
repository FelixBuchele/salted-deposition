// Adapted from the original SDBD source by Nikola Topic (2015).
// See Topic & Pöschel, JCP 308 (2016), DOI: 10.1016/j.jcp.2015.12.052.
#include "matrix.hpp"
#include "../invariant.hpp"
namespace sdbd::detail {
matrix::matrix(const vtype& v1,const vtype& v2,const vtype& v3)//init matrix from three vectors
{
	m11=v1.x();m12=v2.x();m13=v3.x();
	m21=v1.y();m22=v2.y();m23=v3.y();
	m31=v1.z();m32=v2.z();m33=v3.z();
}

vtype ma(const matrix& m, const vtype& v)
{
	vtype res;
	res.x()=m.m11*v.x()+m.m12*v.y()+m.m13*v.z();
	res.y()=m.m21*v.x()+m.m22*v.y()+m.m23*v.z();
	res.z()=m.m31*v.x()+m.m32*v.y()+m.m33*v.z();
	return res;
}
}
