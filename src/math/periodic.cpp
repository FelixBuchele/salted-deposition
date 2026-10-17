// Adapted from the original SDBD source by Nikola Topic (2015).
// See Topic & Pöschel, JCP 308 (2016), DOI: 10.1016/j.jcp.2015.12.052.
#include "periodic.hpp"
#include "../invariant.hpp"
namespace sdbd::detail {
int imod(int i,int n)
{
	if(i>=0) return i%n;
	else return (i%n+n)%n;// Keep exact negative multiples at zero.
}

long double rmod(long double x,long double l)
{
	if(x>=0) return std::fmod(x,l);
	else return std::fmod(std::fmod(x,l)+l,l);
}


long double moddist(long double x, long double l)
{
	return rmod(0.5*l+x,l)-0.5*l;
}

vtype vmod(const vtype& r,long double lx,long double ly)
{
	vtype res;
	res.x()=rmod(r.x(),lx);
	res.y()=rmod(r.y(),ly);
	res.z()=r.z();
	return res;
}

vtype vmod(const vtype& ri,const vtype& rb,long double lx,long double ly)
{
	vtype res;
	res.x()=rb.x()+moddist(ri.x()-rb.x(),lx);
	res.y()=rb.y()+moddist(ri.y()-rb.y(),ly);
	res.z()=ri.z();
	return res;
}
}
