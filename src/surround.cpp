// Adapted from the original SDBD source by Nikola Topic (2015).
// See Topic & Pöschel, JCP 308 (2016), DOI: 10.1016/j.jcp.2015.12.052.
#include "surround.hpp"
#include "engine.hpp"
namespace sdbd::detail {
void surroundtype::init_grid(const Engine& engine)
{
	box.clear();
	box.resize(engine.nx);
	izmax.resize(engine.nx);
	for(int i=0;i<engine.nx;i++) {box[i].resize(engine.ny);izmax[i].resize(engine.ny);}
	for(int i=0;i<engine.nx;i++) 
	for(int j=0;j<engine.ny;j++)
	{
		box[i][j].resize(engine.nz);
		izmax[i][j]=0;
	}
}
}
