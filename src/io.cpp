// Adapted from the original SDBD source by Nikola Topic (2015).
// See Topic & Pöschel, JCP 308 (2016), DOI: 10.1016/j.jcp.2015.12.052.
#include "engine.hpp"
#include "cluster.hpp"
#include "surround.hpp"
#include "geometry_validation.hpp"
namespace sdbd::detail {
namespace {
// Open an output stream with fail/bad exceptions enabled; throw a filename-qualified diagnostic on failure.
void open_output(std::ofstream& stream, const std::string& filename,
                 std::ios::openmode mode = std::ios::out) {
    stream.exceptions(std::ios::failbit | std::ios::badbit);
    try { stream.open(filename, mode); }
    catch (const std::ios_base::failure& error) {
        throw std::runtime_error("Cannot open output " + filename + ": " + error.what());
    }
}
}
void Engine::load_parameters(const char* inputfile) {
    std::ifstream file(inputfile);
    if (!file) throw std::runtime_error(std::string("Cannot open input: ") + inputfile);
    const auto read = [&](const char* label, auto& value) {
        std::string actual;
        if (!(file >> actual >> value) || actual != label)
            throw std::runtime_error(std::string("Invalid input field: ") + label);
    };
    read("N:", N);
    read("nx:", nx); read("ny:", ny); read("nz:", nz);
    read("dx:", dx); read("dy:", dy); read("dz:", dz);
    read("xo:", xo); read("yo:", yo); read("zo:", zo);
    read("phi:", phi);
    read("cmsfile:", cmsfile);
    read("positonsfile:", positionsfile); // Original spelling is part of the format.
    read("detailsfile:", detailsfile);
    std::string extra;
    if (file >> extra) throw std::runtime_error("Unexpected trailing input field: " + extra);
    validate_parameters();
}

void Engine::validate_parameters() {
    geometry_detail::ValidationEnvironment environment;
    if (N <= 0 || nx <= 0 || ny <= 0 || nz <= 0 ||
        !std::isfinite(dx) || !std::isfinite(dy) || !std::isfinite(dz) ||
        dx <= 0 || dy <= 0 || dz <= 0 || !std::isfinite(xo) ||
        !std::isfinite(yo) || !std::isfinite(zo) || zo >= 0 ||
        !std::isfinite(phi) || phi <= 0)
        throw std::invalid_argument("Invalid grid, origin, particle count or phi in input");
    lx=nx*dx; ly=ny*dy;
    if (!std::isfinite(lx) || !std::isfinite(ly) || !std::isfinite(zo+nz*dz) || zo+nz*dz <= 0)
        throw std::invalid_argument("Invalid domain bounds");
    // The legacy periodic cell is [0,lx) x [0,ly); other origins are unsupported.
    if (xo != 0 || yo != 0) throw std::invalid_argument("Horizontal grid origin must be zero");
    if (dx!=dy || dy!=dz || std::floor(dx)!=dx)
        throw std::invalid_argument("Supported grids require equal integer spacings");
    if (static_cast<long double>(nx)*ny*nz>20'000'000)
        throw std::invalid_argument("Grid exceeds supported limit of 20000000 cells");
    pxb=0; pxe=lx; pyb=0; pye=ly;
    geometry_detail::arithmetic_range();
}

vtype Engine::cmass(const std::vector<vtype>& r, const std::vector<long double>& R)
{
	int i;long double mc=0,r3;vtype res=null;
	for(i=0;i<r.size();i++)
	{
		r3=R[i]*R[i]*R[i];
		res=res+r3*(r[i]-r[0]);
		mc=mc+r3;
	}
	res=(1/mc)*res;res=res+r[0];
	return res;
}

int Engine::checkinters(surroundtype& su)
{
	SDBD_REQUIRE(su.particle_ids.size()==su.r.size());
	int count=0;
	for(int ix=0;ix<nx;ix++)
	for(int iy=0;iy<ny;iy++)
	for(int iz=0;iz<nz;iz++)
	for(int ii=0;ii<su.box[ix][iy][iz].size();ii++)
	{
		int i=su.box[ix][iy][iz][ii];
		vtype r=su.r[i];
		long double R=su.R[i];
		int ixl=std::floor((r.x()-R-Rmax-xo)/dx);
		int ixu=std::floor((r.x()+R+Rmax-xo)/dx);

		int iyl=std::floor((r.y()-R-Rmax-yo)/dy);
		int iyu=std::floor((r.y()+R+Rmax-yo)/dy);

		int izl=std::floor((r.z()-R-Rmax-zo)/dz);if(izl<0) izl=0;
		int izu=std::floor((r.z()+R+Rmax-zo)/dz);if(izu>=nz) izu=nz-1;

		for(int ixnt=ixl;ixnt<=ixu;ixnt++)
		for(int iynt=iyl;iynt<=iyu;iynt++)
		for(int iznt=izl;iznt<=izu;iznt++)
		{
			int ixn=imod(ixnt,nx);
			int iyn=imod(iynt,ny);
			int izn=iznt; 
			for(int jj=0;jj<su.box[ixn][iyn][izn].size();jj++)
			{
				int j=su.box[ixn][iyn][izn][jj];
				vtype rj=vmod(su.r[j],r,lx,ly);
				long double Rj=su.R[j];
				if(i!=j && su.particle_ids[i]!=su.particle_ids[j]) if(norm(r-rj)-R-Rj<-60*1e-8) {ds(i,j,norm(r-rj)-R-Rj);throw std::runtime_error("Intersection in deposited spheres");}
			}
		}
	}
	return count/2;
}


bool Engine::checkinters(const clustertype& c,const surroundtype& s)
{
	SDBD_REQUIRE(s.particle_ids.size()==s.r.size());
	for(int i=0;i<c.r.size();i++)
	{
		vtype rp=c.r[i];
		long double R=c.R[i];
		int ixl=std::floor((rp.x()-(R+Rmax)-xo)/dx);
		int ixu=std::floor((rp.x()+(R+Rmax)-xo)/dx);

		int iyl=std::floor((rp.y()-(R+Rmax)-yo)/dy);
		int iyu=std::floor((rp.y()+(R+Rmax)-yo)/dy);

		int izl=std::floor((rp.z()-(R+Rmax)-zo)/dz);if(izl<0) izl=0;
		int izu=std::floor((rp.z()+(R+Rmax)-zo)/dz);if(izu>=nz) izu=nz-1;

		for(int ixnt=ixl;ixnt<=ixu;ixnt++)
		for(int iynt=iyl;iynt<=iyu;iynt++)
		for(int iznt=izl;iznt<=izu;iznt++)
		{
			int ixn=imod(ixnt,nx);
			int iyn=imod(iynt,ny);
			int izn=iznt; 
			for(int jj=0;jj<s.box[ixn][iyn][izn].size();jj++)
			{
				int j=s.box[ixn][iyn][izn][jj];
				{
					vtype rj=vmod(s.r[j],rp,lx,ly);
					if (c.particle_id!=s.particle_ids[j] && norm(c.r[i]-rj)<c.R[i]+s.R[j]-60.*1e-6) {ds("Intersection: ",i,j,norm(c.r[i]-rj)-c.R[i]-s.R[j]);return true;}
				}
			}
		}
	}
	return false;
}


void Engine::find_cluster_contacts(long double eps_contact,const clustertype& c,const surroundtype& su,std::vector<vtype>& rps,std::vector<vtype>& rds,std::vector<int>& p,std::vector<int>& d)
{
	for(int i=0;i<c.r.size();i++)
	{
        if (!c.collision_active(i)) continue;
		long double lcd=c.R[i]+Rmax;
		vtype r=c.r[i];
		int ixl=std::floor((r.x()-lcd-xo)/dx);
		int ixu=std::floor((r.x()+lcd-xo)/dx);

		int iyl=std::floor((r.y()-lcd-yo)/dy);
		int iyu=std::floor((r.y()+lcd-yo)/dy);

		int izl=std::floor((r.z()-lcd-zo)/dz);if(izl<0) izl=0;
		int izu=std::floor((r.z()+lcd-zo)/dz);if(izu>=nz) izu=nz-1;

		for(int ixnt=ixl;ixnt<=ixu;ixnt++)
		for(int iynt=iyl;iynt<=iyu;iynt++)
		for(int iznt=izl;iznt<=izu;iznt++)
		{
			int ixn=imod(ixnt,nx);
			int iyn=imod(iynt,ny);
			int izn=iznt; 
			for(int jj=0;jj<su.box[ixn][iyn][izn].size();jj++)
			{
				int j=su.box[ixn][iyn][izn][jj];
                if (!su.collision_active[j]) continue;
				{
					vtype rj=vmod(su.r[j],r,lx,ly);
					if (norm(r-rj)<c.R[i]+su.R[j]+eps_contact) 
					{
                        // Wrapped grid bounds can visit the same cell more than
                        // once. One sphere pair is still one physical contact.
                        bool recorded=false;
                        for (std::size_t contact=0;contact<p.size();++contact)
                            if (p[contact]==j && d[contact]==i) { recorded=true;break; }
                        if (recorded) continue;
						rps.push_back(rj);	
						rds.push_back(r);	
						p.push_back(j);
						d.push_back(i);
					}
				}
			}
		}
	}
}


void Engine::unload(clustertype& c,surroundtype& su)
{
	su.particle_centers.push_back(c.cm);
    if (c.orientation) su.particle_orientations.push_back(*c.orientation);
	int susize=su.r.size();
	for(int i=0;i<c.r.size();i++) 
	{
		int ii=susize+i;
        unloadsphere(vmod(c.r[i],lx,ly),ii,c.R[i],su,c.particle_id,c.collision_active(i));
	}
}

void Engine::unloadsphere(const vtype& v,int ii,long double R,surroundtype& su,std::size_t particle_id,bool collision_active)
{
	su.r.push_back(v);
	su.R.push_back(R);
    su.particle_ids.push_back(particle_id);
    su.collision_active.push_back(collision_active);

	int ix=std::floor(v.x()-xo)/dx;
	int iy=std::floor(v.y()-yo)/dy;
	int iz=std::floor(v.z()-zo)/dz;
	if (!(0<=ix&&ix<nx&&0<=iy&&iy<ny&&0<=iz&&iz<nz))
        throw std::out_of_range("Deposited sphere exceeds grid bounds");

	su.box[ix][iy][iz].push_back(ii);
	if(iz>su.izmax[ix][iy]) su.izmax[ix][iy]=iz;
}

void Engine::output(surroundtype& su,std::string fname) 
{
	int i;std::ofstream myfile;
	open_output(myfile, fname);
	myfile.precision(15);
	for(i=0;i<su.r.size();i++) {myfile<<su.r[i].x()<<" "<<su.r[i].y()<<" "<<su.r[i].z()<<" "<<su.R[i]<<'\n';}
	myfile.close();
}


void Engine::writedetails(int size,int intersect,std::string name)
{
	std::ofstream myfile;
	
	open_output(myfile, name, std::ios_base::app);
	myfile<<"  Number of complex particles:  "<<size<<"  Number of intersections  "<<intersect<<'\n';
	myfile.close();
}


void Engine::write2cm(long double z,std::string name)
{
	std::ofstream myfile;
	open_output(myfile, name, std::ios::app);
	myfile.precision(15);
	myfile<<2*z<<'\n';
	myfile.close();
}
}
