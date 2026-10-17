// Adapted from the original SDBD source by Nikola Topic (2015).
// See Topic & Pöschel, JCP 308 (2016), DOI: 10.1016/j.jcp.2015.12.052.
#include "rotation.hpp"
#include "../invariant.hpp"
namespace sdbd::detail {
rotinftype rotn(const rotinftype& rotinf1,const rotinftype& rotinf2)
{
	SDBD_REQUIRE(rotinf1.n.x()==rotinf2.n.x()&&rotinf1.n.y()==rotinf2.n.y()&&rotinf1.n.z()==rotinf2.n.z());
	rotinftype res;
	res.cc=rotinf1.cc*rotinf2.cc-rotinf1.ss*rotinf2.ss;	
	res.ss=rotinf1.cc*rotinf2.ss+rotinf1.ss*rotinf2.cc;
	res.n=rotinf1.n;

	return res;
}
}
