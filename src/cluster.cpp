// Adapted from the original SDBD source by Nikola Topic (2015).
// See Topic & Pöschel, JCP 308 (2016), DOI: 10.1016/j.jcp.2015.12.052.
#include "cluster.hpp"
#include "engine.hpp"
#include <array>
namespace sdbd::detail {


vtype clustertype::cmass()
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




clustertype::clustertype(Engine& engine, const ParticleGeometry& geometry, std::size_t id)
    : particle_id(id), explicit_center_(geometry.has_explicit_center()), engine_(&engine) {
    r.reserve(geometry.spheres().size());
    R.reserve(geometry.spheres().size());
    for (const auto& sphere : geometry.spheres()) {
        r.push_back(sphere.center);
        R.push_back(sphere.radius);
    }
    cm=geometry.center_of_mass();
    if (!engine.stlposefile.empty()) orientation=geometry.stl_source()->orientation;
    collision_active_.assign(r.size(),false);
    for (const auto index : geometry.collision_spheres()) collision_active_[index]=true;
}

void clustertype::shift_horizontal(long double x,long double y,bool legacy_recompute_center) {
    if (legacy_recompute_center && explicit_center_) throw std::invalid_argument("Cannot overwrite supplied physical COM with sphere weighting");
    for (auto& point : r) { point.x()+=x; point.y()+=y; }
    if (legacy_recompute_center) cm=cmass();
    else cm+=vtype(x,y,0);
}

bool clustertype::settle(surroundtype& su)
{
	bool collapsed=false;
	int countev=1;
	vtype cmev(cm);

	while(!stable)
	{

		countev++;
		vtype cmp=cm;

		switch(contacts)//propagates cluster by one event
		{
			case 0: fall(su);break;//propagate when cluster is falling
			case 1: rotate1(su);break;//propagate when cluster is rotating on one contact
			case 2: rotate2(su);break;//propagate when cluster is rotating on two contacts
		}

		long double dcm1=norm(cm-cmp);

		long double eps_collapse=60*1e-10;
		long double eps_contact=60*1e-8;
		if(dcm1<eps_collapse&&engine_->particle!=0)//for spheres no need to detect collapse 
		{
			//collapse detected, going through collapse

			std::vector<vtype> rps;
			std::vector<vtype> rds;
			std::vector<int> p;
			std::vector<int> d;
			engine_->find_cluster_contacts(eps_contact,*this,su,rps,rds,p,d);//finds contacts by proximity

			// Repeated fixed-sphere supports do not define distinct rotation axes.
            // Preserve the supplied implementation's collapse/redrop criterion.
			bool duplicate_exists=false;
			for(std::size_t i=0;i+1<p.size();i++)
			for(int j=i+1;j<p.size();j++)
			if(p[i]==p[j]) duplicate_exists=true;

			if(!duplicate_exists)//if there are no duplicates find steepest descent trajectory
			{
				rotattempt att;
				find_steepest_descent(rps,rds,p,d,cm,att);//find steepest descent with closed contacts considered
				//update state
				contacts=att.contacts;
				stable=att.stable;
				p1=att.p1;
				d1=att.d1;
				p2=att.p2;
				d2=att.d2;
				p3=att.p3;
				d3=att.d3;
				td=att.td;
				tp=att.tp;
				f=att.f;
				//particle stays where it is so no need to change cm and r[i]
			}
			else
			{
				// Repeated fixed support: leave retry/freeze policy to the caller.
				stable=true;
				collapsed=true;//considered as collapsed, in order to be redropped, if chosen
			}

		}

		if(countev%1000000==0)//parabolic collapse fixing
		{
			if(norm(cmev-cm)<60) 
			{
				stable=true;
				collapsed=true;
				ds("Moved less than a particle diameter in 1000000 events");
			}
			else cmev=cm;
		}

		/*****************  Important: periodic boundary conditions treatment  *****************/
		periodic(engine_->lx,engine_->ly);
	}
	return collapsed;
}

void clustertype::periodic(long double period_x,long double period_y)
{
	vtype v;
	v.z()=0;

	if(cm.x()<engine_->pxb) v.x()=period_x;
	if(cm.x()>engine_->pxb&&cm.x()<engine_->pxe) v.x()=0;
	if(cm.x()>engine_->pxe) v.x()=-period_x;

	if(cm.y()<engine_->pyb) v.y()=period_y;
	if(cm.y()>engine_->pyb&&cm.y()<engine_->pye) v.y()=0;
	if(cm.y()>engine_->pye) v.y()=-period_y;

	movecluster(v);
}

void clustertype::movecluster(const vtype& v)
{
	for(int i=0;i<r.size();i++) r[i]=r[i]+v;
	cm=cm+v;
}


void clustertype::fall(surroundtype& su)
{
	std::vector<fallattempt> atvec;fallattempt att;

	bfall(su,atvec);//add events for collisions with fixed spheres
	gfall(atvec);//add events for collisions with the ground

	if(atvec.size()==0) {throw std::runtime_error("fall: no events");}//must be at least one future event
	att=*std::max_element(atvec.begin(),atvec.end());//choose the highest future event

	step0(att,su);//move the cluster
}

void clustertype::bfall(surroundtype& su,std::vector<fallattempt>& atvec)
{
	int j,k;long double dlt;fallattempt att;
	bool contact;long double zfc;
	vtype rp;
	for(k=0;k<r.size();k++)//looop through all cluster particles
	{
        if (!collision_active(k)) continue;
		zfc=engine_->h0;contact=false;
		int ixl,ixr,iyl,iyr,izu;engine_->findindex(r[k],R[k],ixl,ixr,iyl,iyr,izu);
		for(int it=ixl;it<=ixr;it++)//loop through boxes containing candidates
		for(int nt=iyl;nt<=iyr;nt++)
		{
			int i=imod(it,engine_->nx);
			int n=imod(nt,engine_->ny);
			for(int m=su.izmax[i][n];m>=0&&(!contact||(contact&&m>=std::floor((zfc-R[k]-engine_->Rmax-engine_->zo)/engine_->dz)));m--)
			for(int s=0;s<su.box[i][n][m].size();s++)
			{
				j=su.box[i][n][m][s];
                if (!su.collision_active[j]) continue;
				rp=vmod(su.r[j],cm,engine_->lx,engine_->ly);
				if(k!=td||j!=tp)//to avoid recontact when conacts was lost by rolling to horizontal
				if(engine_->pcontact0(r[k],R[k],rp,su.R[j],dlt)) //test for collision by translation 
				if(r[k].z()-dlt>rp.z())//gradient check
				{
					if(dlt>0&&!contact) {contact=true;zfc=r[k].z()-dlt;}
					if(dlt>0&&r[k].z()-dlt>zfc) zfc=r[k].z()-dlt;

					att.contacts=1;
					att.stable=false;
					att.p1=j;
					att.d1=k;
					att.p2=-1;
					att.d2=-1;
					att.p3=-1;
					att.d3=-1;
					att.tp=-1;
					att.td=-1;
					att.f=0;
					att.cm=cm;
					att.cm.z()=att.cm.z()-dlt;
					if(att.cm.z()<cm.z()+engine_->epsz) {atvec.push_back(att);}
				}
			}
		}
	}
}

void clustertype::gfall(std::vector<fallattempt>& atvec)
{
	int k;long double dlt;fallattempt att;
	for(k=0;k<r.size();k++) 
	{
        if (!collision_active(k)) continue;
		dlt=r[k].z()-R[k];
		att.stable=true;
		att.contacts=1;
		att.d1=k;
		att.p1=-2;
		att.d2=-1;
		att.p2=-1;
		att.d3=-1;
		att.p3=-1;
		att.td=-1;
		att.tp=-1;
		att.f=0;
		att.cm=cm;
		att.cm.z()=att.cm.z()-dlt;
		if(att.cm.z()<cm.z()+engine_->epsz) atvec.push_back(att);
	}
}

void clustertype::step0(fallattempt att,surroundtype& su)
{
	d1=att.d1;
	p1=att.p1;
	d2=att.d2;
	p2=att.p2;
	d3=att.d3;
	p3=att.p3;
	td=att.td;
	tp=att.tp;
	f=att.f;
	contacts=att.contacts;
	stable=att.stable;
	translate(cm.z()-att.cm.z(),d1,p1,su);
}

void clustertype::translate(long double dlt,int d1,int p1,surroundtype& su)
{
	int i;
	for(i=0;i<r.size();i++) 
	{
		r[i].z()-=dlt;
	}
	cm.z()-=dlt;

	if(p1>=0&&d1>=0)
	{
		//correction due to loss of precision due to use of dlt
		vtype rp1=su.r[p1];
		vtype rd1l=vmod(r[d1],rp1,engine_->lx,engine_->ly);

		long double dl=su.R[p1]+R[d1];
		long double ds2=std::pow(rd1l.x()-rp1.x(),2)+std::pow(rd1l.y()-rp1.y(),2);
		long double underroot=1-ds2/(dl*dl);
		if(underroot>engine_->eps) 
		{
			long double ovrl=norm(rp1-rd1l)-dl;
			r[d1]=r[d1]-ovrl*1/std::sqrt(underroot)*ez;
			rd1l=vmod(r[d1],rp1,engine_->lx,engine_->ly);
		}
	}
}

void clustertype::rotate1(surroundtype& su)
{
	std::vector<rotattempt> atvec;rotattempt att;
	vtype rp_temp=vmod(su.r[p1],cm,engine_->lx,engine_->ly);

	//correction of the case where the center of mass is above the center of fixed particle
	bool cm_on_top_of_fixed_particle=std::fabs(cm.x()-rp_temp.x())<60*1e-6&&std::fabs(cm.y()-rp_temp.y())<60*1e-6;
	vtype sh;

	if(cm_on_top_of_fixed_particle)//moves p1 particle slightly
	{	
        // Two uniform horizontal draws reproduce the supplied code's square
        // perturbation; its contradictory rejection predicate never repeated.
        sh.z()=0;
        sh.x()=engine_->ranfs();
        sh.y()=engine_->ranfs();
		sh=engine_->eps*sh;
		su.r[p1]=su.r[p1]+sh;
	}
	int p1old=p1;

	vtype rp=vmod(su.r[p1],cm,engine_->lx,engine_->ly);

	steeperdphi(rp,atvec);//dphi event
	brot1(su,atvec);//add events for collisions with fixed spheres
	hrot(rp,atvec);//add event for rotation to horizontal 
	bottomrot1(rp,atvec);//add event for rotation to hanging
	grot1(rp,atvec);//add event for contact with the ground
	if(atvec.size()==0) {throw std::runtime_error("rotate1: no events");}

	att=*std::max_element(atvec.begin(),atvec.end());
	step1(att,rp);//moves cluster

	if(cm_on_top_of_fixed_particle) su.r[p1old]=su.r[p1old]-sh;
}

void clustertype::steeperdphi(const vtype& rp,std::vector<rotattempt>& atvec)
{
	if(r[d1].z()<rp.z())//if cluster particle is below fixed one then cluster is frustrated situation
	{		
		rotinftype rotinf=engine_->rotinfphi;
		vtype cmr=cm-rp;
		vtype n=v(ez,cmr);n=(1/norm(n))*n;
		rotinf.n=n;
		
		vtype cmc=rot(cm,rp,rotinf.cc,rotinf.ss,rotinf.n);//rotated center of mass
		vtype rdr=rot(r[d1],rp,rotinf.cc,rotinf.ss,rotinf.n);//rotated cluster particle d1

		if(engine_->rightside(cm-rp,cmc-rp,rotinf.n))//checks whether the cluster is on the same side of the plane with normal ez x n, please see paper
		{
			if(rdr.z()<rp.z())//check whether after rotation by dphi angle cluster sphere is still below fixed sphere
			{
				rotattempt att;
				att.contacts=0;
				att.stable=false;
				att.d1=-1;
				att.p1=-1;
				att.d2=-1;
				att.p2=-1;
				att.d3=-1;
				att.p3=-1;
				att.td=-1;
				att.tp=-1;
				att.f=0;
				att.cm=cmc;
				att.rotinf=rotinf;
				atvec.push_back(att);
				SDBD_REQUIRE(cmc.z()<cm.z()+engine_->epsz);
			}
		}
	}
}


void clustertype::brot1(surroundtype& su,std::vector<rotattempt>& atvec)
{
	rotattempt att;
	int j;
	vtype rp=vmod(su.r[p1],cm,engine_->lx,engine_->ly);

//find
	long double rmax=0;
	long double lmax=0;
	for(int k=0;k<r.size();k++)//find the rotation angle between test positions of the cluster such that no events are missed
	{
        if (!collision_active(k)) continue;
		long double rt=engine_->distance_from_axis(r[k],rp,cm);
		if(rt>rmax) 
		{
			rmax=rt;
			lmax=2*(engine_->Rmax+R[k]);
		}
	}
	long double dephi;

	if(rmax<lmax) dephi=engine_->l_pi/3.;
	else dephi=2*std::asin(lmax/(2*rmax));//angle increment

	long double phic=0;
	long double phidet=engine_->l_pi;
	do//do while loop until first contact has been detected without doubt
	{
		for(int k=0;k<r.size();k++)//loop through cluster particles
		{
        if (!collision_active(k)) continue;
			vtype point=engine_->rotpoint(r[k],rp,cm,phic);//position of the cluster particle according to current test angle
			long double axisr=engine_->distance_from_axis(r[k],rp,cm);
			long double l=2*axisr*std::sin(dephi/2);
			long double max_cillinder=R[k]+engine_->Rmax+axisr-std::sqrt(axisr*axisr-l*l/4.);
			long double radius=std::sqrt(l*l/4.+max_cillinder*max_cillinder);//radius of the search volume

			vtype rref=point;//center of spheres giving the neighbours

			//testing collision with candidate spheres 
			int ixl=std::floor((rref.x()-radius-engine_->xo)/engine_->dx);
			int ixu=std::floor((rref.x()+radius-engine_->xo)/engine_->dx);

			int iyl=std::floor((rref.y()-radius-engine_->yo)/engine_->dy);
			int iyu=std::floor((rref.y()+radius-engine_->yo)/engine_->dy);

			int izl=std::floor((rref.z()-radius-engine_->zo)/engine_->dz);if(izl<0) izl=0;
			int izu=std::floor((rref.z()+radius-engine_->zo)/engine_->dz);if(izu>=engine_->nz) izu=engine_->nz-1;

			for(int ixnt=ixl;ixnt<=ixu;ixnt++)
			for(int iynt=iyl;iynt<=iyu;iynt++)
			for(int iznt=izl;iznt<=izu;iznt++)
			{
				int ixn=imod(ixnt,engine_->nx);
				int iyn=imod(iynt,engine_->ny);
				int izn=iznt; 

				for(int jj=0;jj<su.box[ixn][iyn][izn].size();jj++)
				{
					j=su.box[ixn][iyn][izn][jj];
                    if (!su.collision_active[j]) continue; // index of the candidate sphere
					vtype rp2=vmod(su.r[j],cm,engine_->lx,engine_->ly);//periodic conditions applied to candidate sphere
					{
						rotinftype rotinf1,rotinf2,rotinf;//information about rotation
						if(j!=p1)//no need to test for contact 1
						if(!f||k!=td||j!=tp)//to avoid reconctact after grazing loss of contact
						if(engine_->pcontact1(rp,cm,r[k],R[k],rp2,su.R[j],rotinf1,rotinf2))//test for collision
						{
							const std::array rotv{rotinf1,rotinf2};


							for(int ir=0;ir<2;ir++)//consider both solutions for collision
							{
								rotinf=rotv[ir];

								//rotated relevant quantities
								vtype cmc=rot(cm,rp,rotinf.cc,rotinf.ss,rotinf.n);
								vtype rdc=rot(r[d1],rp,rotinf.cc,rotinf.ss,rotinf.n);
								vtype rcc=rot(r[k],rp,rotinf.cc,rotinf.ss,rotinf.n);
								if(engine_->rightside(cm-rp,cmc-rp,rotinf.n)&&engine_->g(rp,cmc,rcc,rp2)<0&&cmc.z()<cm.z()+engine_->epsz)//final check if collision is a valid future event
								{
									std::vector<vtype> rps;
									rps.push_back(rp);rps.push_back(rp2);
									std::vector<vtype> rds;
									rds.push_back(rdc);rds.push_back(rcc);
									std::vector<int> p;
									p.push_back(p1);p.push_back(j);
									std::vector<int> d;
									d.push_back(d1);d.push_back(k);
									find_steepest_descent(rps,rds,p,d,cmc,att);//consider all contacts and find steepest descent trajectory
									att.rotinf=rotinf; 
									atvec.push_back(att);//add event

									//if by numerical error |cos|>1
									if(att.rotinf.cc<-1) att.rotinf.cc=att.rotinf.cc+1e-12;
									if(att.rotinf.cc>1) att.rotinf.cc=att.rotinf.cc-1e-12;
									long double phidetc=std::acos(att.rotinf.cc);//rotation angle for the current collision
									if(phidetc<phidet) phidet=phidetc;//update smallest colision angle
								}//right side if
							}//+- solution for loop
						}
					}
				}
			}		
		}
		phic+=dephi;//increase the angle for the cluster test position
	}
	while((phic-dephi)<phidet);//checks whether the current test angle is too large compared to current collision angle, for earlier collision to be possible  
}

void clustertype::hrot(const vtype& rp,std::vector<rotattempt>& atvec)
{
	rotattempt att;
	rotinftype rotinf;
	engine_->tohorizontal(cm,r[d1],rp,rotinf);
	att.contacts=0;
	att.stable=false;
	att.p1=-1;
	att.d1=-1;
	att.p2=-1;
	att.d2=-1;
	att.p3=-1;
	att.d3=-1;
	att.td=d1;
	att.tp=p1;
	att.f=0;

	att.cm=rot(cm,rp,rotinf.cc,rotinf.ss,rotinf.n);
	att.rotinf=rotinf;

	vtype cmc=rot(cm,rp,rotinf.cc,rotinf.ss,rotinf.n);
	if(engine_->rightside(cm-rp,cmc-rp,rotinf.n)&&att.cm.z()<cm.z()+engine_->epsz) atvec.push_back(att);
}

void clustertype::grot1(const vtype& rp,std::vector<rotattempt>& atvec)
{
	rotattempt att;
	rotinftype rotinf1,rotinf2;
	for(int i=0;i<r.size();i++)
	{
        if (!collision_active(i)) continue;
		if(engine_->toground1(rp,cm,r[i],R[i],rotinf1,rotinf2))//finds position (if any) where cluster sphere i makes contact with the ground
		{
			vtype cmc1=rot(cm,rp,rotinf1.cc,rotinf1.ss,rotinf1.n);
			vtype cmc2=rot(cm,rp,rotinf2.cc,rotinf2.ss,rotinf2.n);
			if(engine_->rightside(cm-rp,cmc1-rp,rotinf1.n))
			{
				att.contacts=2;
				att.stable=true;
				att.d1=d1;
				att.p1=p1;
				att.d2=i;
				att.p2=-2;
				att.d3=-1;
				att.p3=-1;
				att.td=-1;
				att.tp=-1;
				att.f=0;
				att.cm=cmc1;
				att.rotinf=rotinf1;
				if(att.cm.z()<cm.z()+engine_->epsz) atvec.push_back(att);
			}

			if(engine_->rightside(cm-rp,cmc2-rp,rotinf2.n))
			{
				att.contacts=2;
				att.stable=true;
				att.d1=d1;
				att.p1=p1;
				att.d2=i;
				att.p2=-2;
				att.d3=-1;
				att.p3=-1;
				att.td=-1;
				att.tp=-1;
				att.f=0;
				att.cm=cmc2;
				att.rotinf=rotinf2;
				if(att.cm.z()<cm.z()+engine_->epsz) atvec.push_back(att);
			}
		}
	}
}

void clustertype::bottomrot1(const vtype& rp,std::vector<rotattempt>& atvec)
{
	rotattempt att;
	rotinftype rotinf;
	engine_->tobelow1(cm,rp,rotinf);//finds amount of rotation such that the center of mass is below the fixed particle

	att.contacts=30;
	att.stable=true;
	att.p1=-1;
	att.d1=-1;
	att.p2=-1;
	att.d2=-1;
	att.p3=-1;
	att.d3=-1;
	att.cm=rot(cm,rp,rotinf.cc,rotinf.ss,rotinf.n);
	att.rotinf=rotinf;
	atvec.push_back(att);
	SDBD_REQUIRE(att.cm.z()<=cm.z()+engine_->epsz);
}

void clustertype::step1(const rotattempt& att,const vtype& rp)
{
	contacts=att.contacts;
	stable=att.stable;
	p1=att.p1;
	d1=att.d1;
	p2=att.p2;
	d2=att.d2;
	p3=att.p3;
	d3=att.d3;
	td=att.td;
	tp=att.tp;
	f=att.f;
	cm=rot(cm,rp,att.rotinf.cc,att.rotinf.ss,att.rotinf.n);
	for(int i=0;i<r.size();i++) r[i]=rot(r[i],rp,att.rotinf.cc,att.rotinf.ss,att.rotinf.n);
    if (orientation) *orientation=Quaternion::from_axis(att.rotinf.n,att.rotinf.cc,att.rotinf.ss)* *orientation;
}

void clustertype::rotate2(surroundtype& su)
{
	std::vector<rotattempt> atvec;rotattempt att;

	vtype rp1_temp=vmod(su.r[p1],cm,engine_->lx,engine_->ly);
	vtype rp2_temp=vmod(su.r[p2],cm,engine_->lx,engine_->ly);

	vtype n=rp2_temp-rp1_temp;n=(1/norm(n))*n;
	vtype cmr=cm-rp1_temp;cmr=(1/norm(cmr))*cmr;
	vtype vcmn=v(n,cmr);vcmn=(1/norm(vcmn))*vcmn;

	vtype sh;
	bool is_middle=std::fabs(s(ez,vcmn))<60*1e-9;
	if(is_middle) 
	{
		long double sh_val;
		do
		{
			sh_val=engine_->ranfs();
		}
		while(std::fabs(sh_val)<0.5);	
		sh=60*1e-9*sh_val*vcmn;
		su.r[p1]=su.r[p1]+sh;//moves fixed particle by small amount when cm is exactly on top of rotation axis
	}
	int p1old=p1;

	vtype rp1=vmod(su.r[p1],cm,engine_->lx,engine_->ly);
	vtype rp2=vmod(su.r[p2],cm,engine_->lx,engine_->ly);
	brot2(su,atvec);//add events for collisions with fixed spheres
	disconnect(rp1,rp2,atvec,su);//add events for the case Fig6a of the paper
	steeperdphi2(rp1,rp2,atvec);//add event for "frustrated" situation
	hrot2(rp1,rp2,atvec);//add events for the case Fig6b of the paper
	bottomrot2(rp1,rp2,atvec);//add events for the case Fig6c of the paper
	grot2(rp1,rp2,atvec);//add events for contact with the ground

	if(atvec.size()==0) {throw std::runtime_error("rotate2: no events");}

	att=*std::max_element(atvec.begin(),atvec.end());//choose earliest event
	step2(att,su,rp1);//propagate cluster by chosen event

	if(is_middle) su.r[p1old]=su.r[p1old]-sh;
}

void clustertype::steeperdphi2(const vtype& rp1,const vtype& rp2,std::vector<rotattempt>& atvec)
{
	//it is guarantied that rp1 corresponds to p1 and rp2 to p2, see function rotate2()
	std::vector<vtype> rp;rp.push_back(rp1);rp.push_back(rp2);
	std::vector<vtype> rd;rd.push_back(r[d1]);rd.push_back(r[d2]);
	std::vector<int> p;p.push_back(p1);p.push_back(p2);
	std::vector<int> d;d.push_back(d1);d.push_back(d2);

	rotattempt att;
	find_steepest_descent(rp,rd,p,d,cm,att);//find steepest descent trajectory

	if(att.contacts<contacts)//is it possible to move with less contacts?
	{
		[[maybe_unused]] vtype cmr=cm-rp1;vtype n=rp2-rp1;n=(1/norm(n))*n;
		if(s(ez,v(n,cm-rp1))>0) n=-1*n;

		rotinftype rotinf=engine_->rotinfphi;
		rotinf.n=n;

		vtype cmc=rot(cm,rp1,rotinf.cc,rotinf.ss,rotinf.n);//rotate center of mass of the cluster by small angle

		if(engine_->rightside(cm-rp1,cmc-rp1,rotinf.n))
		{
			//vectors rp and rd
			std::vector<vtype> rp;rp.push_back(rp1);rp.push_back(rp2);

			vtype rd1r=rot(r[d1],rp1,rotinf.cc,rotinf.ss,rotinf.n);
			vtype rd2r=rot(r[d2],rp1,rotinf.cc,rotinf.ss,rotinf.n);
			std::vector<vtype> rd;rd.push_back(rd1r);rd.push_back(rd2r);
			
			//indices of contacting particles
			std::vector<int> p;p.push_back(p1);p.push_back(p2);
			std::vector<int> d;d.push_back(d1);d.push_back(d2);

			rotattempt att;
			find_steepest_descent(rp,rd,p,d,cmc,att);//find steepest descent after the cluster has been rotated
			
			if(att.contacts<contacts) 
			{
				att.rotinf=rotinf;
				atvec.push_back(att);
				SDBD_REQUIRE(att.cm.z()<=cm.z());
			}
			else //consistency check that rolling on two has the same contacts
			{
				SDBD_REQUIRE(att.stable==stable);
				SDBD_REQUIRE(att.contacts==contacts);
				SDBD_REQUIRE((att.p1==p1&&att.d1==d1)||(att.p2==p1&&att.d2==d1));
			}
		}
	}
	else
	{
		//check of consistency between search for steepest descent in previous and current step
		SDBD_REQUIRE(att.stable==stable);
		SDBD_REQUIRE(att.contacts==contacts);
		SDBD_REQUIRE((att.p1==p1&&att.d1==d1)||(att.p2==p1&&att.d2==d1));
	}
}

void clustertype::brot2(surroundtype& su,std::vector<rotattempt>& atvec)
{
	rotattempt att;
	rotinftype rotinf1,rotinf2,rotinf;

	vtype rp=vmod(su.r[p1],cm,engine_->lx,engine_->ly);
	vtype rp2=vmod(su.r[p2],cm,engine_->lx,engine_->ly);

	long double rmax=0;
	long double lmax=0;
	for(int k=0;k<r.size();k++)
	{
        if (!collision_active(k)) continue;
		long double rt=engine_->distance_from_axis2(r[k],rp,rp2);
		if(rt>rmax) 
		{
			rmax=rt;
			lmax=2*(engine_->Rmax+R[k]);
		}
	}
	long double dephi;

	if(rmax<lmax) dephi=engine_->l_pi/3.;
	else dephi=2*std::asin(lmax/(2*rmax));


	long double phic=0;
	long double phidet=engine_->l_pi;
	do
	{
		for(int k=0;k<r.size();k++)
		{
        if (!collision_active(k)) continue;
			vtype point=engine_->rotpoint2(r[k],rp,rp2,cm,phic);
			long double axisr=engine_->distance_from_axis2(r[k],rp,rp2);
			long double l=2*axisr*std::sin(dephi/2);
			long double max_cillinder=R[k]+engine_->Rmax+axisr-std::sqrt(axisr*axisr-l*l/4.);
			long double radius=std::sqrt(l*l/4.+max_cillinder*max_cillinder);

			vtype rref=point; //center of spheres giving the neighbours
			int ixl=std::floor((rref.x()-radius-engine_->xo)/engine_->dx);
			int ixu=std::floor((rref.x()+radius-engine_->xo)/engine_->dx);

			int iyl=std::floor((rref.y()-radius-engine_->yo)/engine_->dy);
			int iyu=std::floor((rref.y()+radius-engine_->yo)/engine_->dy);

			int izl=std::floor((rref.z()-radius-engine_->zo)/engine_->dz);if(izl<0) izl=0;
			int izu=std::floor((rref.z()+radius-engine_->zo)/engine_->dz);if(izu>=engine_->nz) izu=engine_->nz-1;

			for(int ixnt=ixl;ixnt<=ixu;ixnt++)
			for(int iynt=iyl;iynt<=iyu;iynt++)
			for(int iznt=izl;iznt<=izu;iznt++)
			{
				int ixn=imod(ixnt,engine_->nx);
				int iyn=imod(iynt,engine_->ny);
				int izn=iznt;

				for(int jj=0;jj<su.box[ixn][iyn][izn].size();jj++)
				{
					int j=su.box[ixn][iyn][izn][jj];
                    if (!su.collision_active[j]) continue;
					vtype rp3=vmod(su.r[j],cm,engine_->lx,engine_->ly);
					if(j!=p1&&j!=p2)
					{
						if(engine_->pcontact2(rp,rp2,cm,r[k],R[k],rp3,su.R[j],rotinf1,rotinf2))
						{
//							vector<gradrotattempt> gradevent;
//							gradevent.clear();

							const std::array rotv{rotinf1,rotinf2};


							for(int ir=0;ir<2;ir++)
							{
								rotinf=rotv[ir];
								vtype cmc=rot(cm,rp,rotinf.cc,rotinf.ss,rotinf.n);
								vtype rcr=rot(r[k],rp,rotinf.cc,rotinf.ss,rotinf.n);
								vtype rd1c=rot(r[d1],rp,rotinf.cc,rotinf.ss,rotinf.n);
								vtype rd2c=rot(r[d2],rp,rotinf.cc,rotinf.ss,rotinf.n);
								if(engine_->rightside(cm-rp,cmc-rp,rotinf.n)&&engine_->g2(rp,rp2,cmc,rcr,rp3)<0&&cmc.z()<cm.z()+engine_->epsz)
								{
									std::vector<vtype> rps;
									rps.push_back(rp);rps.push_back(rp2);rps.push_back(rp3);
									std::vector<vtype> rds;
									rds.push_back(rd1c);rds.push_back(rd2c);rds.push_back(rcr);
									std::vector<int> p;
									p.push_back(p1);p.push_back(p2);p.push_back(j);
									std::vector<int> d;
									d.push_back(d1);d.push_back(d2);d.push_back(k);
									find_steepest_descent(rps,rds,p,d,cmc,att);
									att.rotinf=rotinf;
									atvec.push_back(att);
									if(att.rotinf.cc<-1) att.rotinf.cc=att.rotinf.cc+1e-12;
	if(att.rotinf.cc>1) att.rotinf.cc=att.rotinf.cc-1e-12;
									long double phidetc=std::acos(att.rotinf.cc);
									if(phidetc<phidet) phidet=phidetc;
								}//right side if
							}//soulution for loop
						}
					}
				}
			}
		}
		phic+=dephi;
	}
	while((phic-dephi)<phidet);
}

void clustertype::bottomrot2(const vtype& rp1,const vtype& rp2,std::vector<rotattempt>& atvec)
{
	rotattempt att;
	rotinftype rotinf;
	engine_->tobelow2(cm,rp1,rp2,rotinf);
	vtype cmc=rot(cm,rp1,rotinf.cc,rotinf.ss,rotinf.n);
	att.contacts=20;//cluster is stable so this value is irrelevant. Value 20 signals that cluster hangs
	att.stable=true;
	att.p1=p1;
	att.d1=d1;
	att.p2=p2;
	att.d2=d2;
	att.p3=-1;
	att.d3=-1;
	att.tp=-1;
	att.td=-1;
	att.f=0;
	att.cm=cmc;
	att.rotinf=rotinf;
	atvec.push_back(att);
	SDBD_REQUIRE(att.cm.z()<=cm.z()+engine_->epsz); 
	return;
}

void clustertype::grot2(const vtype rp1,const vtype rp2,std::vector<rotattempt>& atvec)
{
	rotattempt att;
	rotinftype rotinf1,rotinf2;
	for(int i=0;i<r.size();i++)
	{
        if (!collision_active(i)) continue;
		if(engine_->toground2(rp1,rp2,cm,r[i],R[i],rotinf1,rotinf2))
		{
			vtype cmc1=rot(cm,rp1,rotinf1.cc,rotinf1.ss,rotinf1.n);
			vtype cmc2=rot(cm,rp1,rotinf2.cc,rotinf2.ss,rotinf2.n);
			if(engine_->rightside(cm-rp1,cmc1-rp1,rotinf1.n))
			{
				att.contacts=3;
				att.stable=true;
				att.d1=d1;
				att.p1=p1;
				att.d2=d2;
				att.p2=p2;
				att.d3=i;
				att.p3=-2;
				att.td=-1;
				att.tp=-1;
				att.f=0;
				att.cm=cmc1;
				att.rotinf=rotinf1;
				if(att.cm.z()<cm.z()+engine_->epsz) atvec.push_back(att);
			}

			if(engine_->rightside(cm-rp1,cmc2-rp1,rotinf2.n))
			{
				att.contacts=3;
				att.stable=true;
				att.d1=d1;
				att.p1=p1;
				att.d2=d2;
				att.p2=p2;
				att.d3=i;
				att.p3=-2;
				att.td=-1;
				att.tp=-1;
				att.f=0;
				att.cm=cmc2;
				att.rotinf=rotinf2;
				if(att.cm.z()<cm.z()+engine_->epsz) atvec.push_back(att);
			}
		}
	}
}

void clustertype::hrot2(const vtype& rp1,const vtype& rp2,std::vector<rotattempt>& atvec)
{
	rotattempt att;
	rotinftype rotinf1,rotinf2;
	if(r[d1].z()>rp1.z()) 
	{
		if(engine_->tohorizontal2(rp1,rp2,cm,r[d1],rotinf1,rotinf2))
		{

			const std::array rotv{rotinf1,rotinf2};

			for(int i=0;i<rotv.size();i++)
			{
				rotinftype rotinf=rotv[i];
				vtype cmc=rot(cm,rp1,rotinf.cc,rotinf.ss,rotinf.n);
				[[maybe_unused]] vtype rd1r=rot(r[d1],rp1,rotinf.cc,rotinf.ss,rotinf.n);			
				[[maybe_unused]] vtype rd2r=rot(r[d2],rp1,rotinf.cc,rotinf.ss,rotinf.n);
				if(engine_->rightside(cm-rp1,cmc-rp1,rotinf.n)&&rd2r.z()<rp2.z())
				{
					att.contacts=0;
					att.stable=false;
					att.p1=-1;
					att.d1=-1;
					att.p2=-1;
					att.d2=-1;
					att.p3=-1;
					att.d3=-1;
					att.tp=p1;
					att.td=d1;
					att.f=0;
					att.cm=cmc;
					att.rotinf=rotinf;
					if(att.cm.z()<cm.z()+engine_->epsz) {atvec.push_back(att);}			
				}
			}
		}
	}

	if(r[d2].z()>rp2.z()) 
	{
		if(engine_->tohorizontal2(rp2,rp1,cm,r[d2],rotinf1,rotinf2))
		{

			const std::array rotv{rotinf1,rotinf2};

			for(int i=0;i<rotv.size();i++)
			{
				rotinftype rotinf=rotv[i];
				vtype cmc=rot(cm,rp2,rotinf.cc,rotinf.ss,rotinf.n);
				[[maybe_unused]] vtype rd1r=rot(r[d1],rp2,rotinf.cc,rotinf.ss,rotinf.n);			
				[[maybe_unused]] vtype rd2r=rot(r[d2],rp2,rotinf.cc,rotinf.ss,rotinf.n);
				if(engine_->rightside(cm-rp2,cmc-rp2,rotinf.n)&&rd1r.z()<rp1.z())
				{
					att.contacts=0;
					att.stable=false;
					att.p1=-1;
					att.d1=-1;
					att.p2=-1;
					att.d2=-1;
					att.p3=-1;
					att.d3=-1;
					att.tp=p2;
					att.td=d2;
					att.f=0;
					att.cm=cmc;
					att.rotinf=rotinf;
					if(att.cm.z()<cm.z()+engine_->epsz) {atvec.push_back(att);}			
				}
			}
		}
	}
}

void clustertype::disconnect(const vtype& rp1,const vtype& rp2, std::vector<rotattempt>& atvec,surroundtype& /*su*/)
{
	rotattempt att;
	rotinftype rotinf1,rotinf2;
	vtype rd1=r[d1];
	vtype rd2=r[d2];

	//treatment of contact 2 
	bool push12=engine_->g(rp1,cm,rd2,rp2)<0;//is contact 2 closing while rolling just on contact 1?
	if(push12&&engine_->pdisconnect(cm,rp1,rp2,rd1,rd2,rotinf1,rotinf2))//if so find position of the cluster where d2 moves tangentially to p2, while rotation is on p1 (see paper Eq. 25)
	{
		//starting solution 1
		vtype cmc1=rot(cm,rp1,rotinf1.cc,rotinf1.ss,rotinf1.n);
		if(engine_->rightside(cm-rp1,cmc1-rp1,rotinf1.n))
		{
			if(positivecircle(rp1,rp2,rd2,rotinf1))//checking the condition (31) of the paper
			{
				//condition (31) not satisfied. Standard treatment.
				att.contacts=1;
				att.stable=false;
				att.p1=p1;
				att.d1=d1;
				att.p2=-1;
				att.d2=-1;
				att.p3=-1;
				att.d3=-1;
				att.tp=p2;
				att.td=d2;
				att.f=1;
				att.cm=cmc1;
				att.rotinf=rotinf1;
				if(att.cm.z()<cm.z()+engine_->epsz) atvec.push_back(att);
			}
			else
			{
				//condition (31) satisfied. Cluster needs to be rotated by dphi before contact can be lost
				engine_->rotinfphi.n=rotinf1.n;
				rotinftype rotinfdphi=rotn(rotinf1,engine_->rotinfphi);
				vtype cmcdph=rot(cm,rp1,rotinfdphi.cc,rotinfdphi.ss,rotinfdphi.n);
				vtype rd2r=rot(rd2,rp1,rotinfdphi.cc,rotinfdphi.ss,rotinfdphi.n);
				if(engine_->rightside(cm-rp1,cmcdph-rp1,rotinfdphi.n)&&engine_->g(rp1,cmcdph,rd2r,rp2)>engine_->epstang)
				{
					att.contacts=1;
					att.stable=false;
					att.p1=p1;
					att.d1=d1;
					att.p2=-1;
					att.d2=-1;
					att.p3=-1;
					att.d3=-1;
					att.tp=-1;
					att.td=-1;
					att.f=1;
					att.cm=cmcdph;
					att.rotinf=rotinfdphi;
					if(att.cm.z()<cm.z()+engine_->epsz) {atvec.push_back(att);}
				}
			}
		}

		//starting solution 2
		vtype cmc2=rot(cm,rp1,rotinf2.cc,rotinf2.ss,rotinf2.n);
		if(engine_->rightside(cm-rp1,cmc2-rp1,rotinf2.n))
		{
			if(positivecircle(rp1,rp2,rd2,rotinf2))
			{
				att.contacts=1;
				att.stable=false;
				att.p1=p1;
				att.d1=d1;
				att.p2=-1;
				att.d2=-1;
				att.p3=-1;
				att.d3=-1;
				att.tp=p2;
				att.td=d2;
				att.f=1;
				att.cm=cmc2;
				att.rotinf=rotinf2;
				if(att.cm.z()<cm.z()+engine_->epsz) atvec.push_back(att);
			}
			else
			{
				engine_->rotinfphi.n=rotinf2.n;
				rotinftype rotinfdphi=rotn(rotinf2,engine_->rotinfphi);
				vtype cmcdph=rot(cm,rp1,rotinfdphi.cc,rotinfdphi.ss,rotinfdphi.n);
				vtype rd2r=rot(rd2,rp1,rotinfdphi.cc,rotinfdphi.ss,rotinfdphi.n);
				if(engine_->rightside(cm-rp1,cmcdph-rp1,rotinfdphi.n)&&engine_->g(rp1,cmcdph,rd2r,rp2)>engine_->epstang)
				{
					att.contacts=1;
					att.stable=false;
					att.p1=p1;
					att.d1=d1;
					att.p2=-1;
					att.d2=-1;
					att.p3=-1;
					att.d3=-1;
					att.tp=-1;
					att.td=-1;
					att.f=1;
					att.cm=cmcdph;
					att.rotinf=rotinfdphi;
					if(att.cm.z()<cm.z()+engine_->epsz) {atvec.push_back(att);}
				}	
			}
		}
	}

	//end of treatment of contact2

	//treatment of contact 1
	bool push21=engine_->g(rp2,cm,rd1,rp1)<0;
	if(push21&&engine_->pdisconnect(cm,rp2,rp1,rd2,rd1,rotinf1,rotinf2))
	{
		//solution 1
		vtype cmc1=rot(cm,rp1,rotinf1.cc,rotinf1.ss,rotinf1.n);
		if(engine_->rightside(cm-rp1,cmc1-rp1,rotinf1.n))
		{
			if(positivecircle(rp2,rp1,rd1,rotinf1))
			{
				att.contacts=1;
				att.stable=false;
				att.p1=p2;
				att.d1=d2;
				att.p2=-1;
				att.d2=-1;
				att.p3=-1;
				att.d3=-1;
				att.tp=p1;
				att.td=d1;
				att.f=1;
				att.cm=cmc1;
				att.rotinf=rotinf1;
				if(att.cm.z()<cm.z()+engine_->epsz) atvec.push_back(att);
			}
			else
			{
				engine_->rotinfphi.n=rotinf1.n;
				rotinftype rotinfdphi=rotn(rotinf1,engine_->rotinfphi);
				vtype cmcdph=rot(cm,rp1,rotinfdphi.cc,rotinfdphi.ss,rotinfdphi.n);
				vtype rd1r=rot(rd1,rp1,rotinfdphi.cc,rotinfdphi.ss,rotinfdphi.n);
				if(engine_->rightside(cm-rp1,cmcdph-rp1,rotinfdphi.n)&&engine_->g(rp2,cmcdph,rd1r,rp1)>engine_->epstang)
				{
					att.contacts=1;
					att.stable=false;
					att.p1=p2;
					att.d1=d2;
					att.p2=-1;
					att.d2=-1;
					att.p3=-1;
					att.d3=-1;
					att.tp=-1;
					att.td=-1;
					att.f=0;
					att.cm=cmcdph;
					att.rotinf=rotinfdphi;
					if(att.cm.z()<cm.z()+engine_->epsz) {atvec.push_back(att);}			
				}
			}
		}

		//solution 2
		vtype cmc2=rot(cm,rp1,rotinf2.cc,rotinf2.ss,rotinf2.n);
		[[maybe_unused]] vtype rd2c2=rot(rd2,rp1,rotinf2.cc,rotinf2.ss,rotinf2.n);
		if(engine_->rightside(cm-rp1,cmc2-rp1,rotinf2.n))
		{
			if(positivecircle(rp2,rp1,rd1,rotinf2))
			{
				att.contacts=1;
				att.stable=false;
				att.p1=p2;
				att.d1=d2;
				att.p2=-1;
				att.d2=-1;
				att.p3=-1;
				att.d3=-1;
				att.tp=p1;
				att.td=d1;
				att.f=1;
				att.cm=cmc2;
				att.rotinf=rotinf2;
				if(att.cm.z()<cm.z()+engine_->epsz) atvec.push_back(att);
			}
			else
			{
				engine_->rotinfphi.n=rotinf2.n;
				rotinftype rotinfdphi=rotn(rotinf2,engine_->rotinfphi);
				vtype cmcdph=rot(cm,rp1,rotinfdphi.cc,rotinfdphi.ss,rotinfdphi.n);
				vtype rd1r=rot(rd1,rp1,rotinfdphi.cc,rotinfdphi.ss,rotinfdphi.n);

				[[maybe_unused]] vtype cmcn=rot(cm,rp1,rotinf2.cc,rotinf2.ss,rotinf2.n);
				[[maybe_unused]] vtype cmcnd=rot(cm,rp1,rotinfdphi.cc,rotinfdphi.ss,rotinfdphi.n);

				if(engine_->rightside(cm-rp1,cmcdph-rp1,rotinfdphi.n)&&engine_->g(rp2,cmcdph,rd1r,rp1)>engine_->epstang)
				{
					att.contacts=1;
					att.stable=false;
					att.p1=p2;
					att.d1=d2;
					att.p2=-1;
					att.d2=-1;
					att.p3=-1;
					att.d3=-1;
					att.tp=-1;
					att.td=-1;
					att.f=0;
					att.cm=cmcdph;
					att.rotinf=rotinfdphi;
					if(att.cm.z()<cm.z()+engine_->epsz) {atvec.push_back(att);}			
				}
			}
		}
	}
}

bool clustertype::positivecircle(const vtype& rp1,const vtype& rp2,const vtype& rd2,const rotinftype& rotinf)
{
	vtype cmr=rot(cm,rp1,rotinf.cc,rotinf.ss,rotinf.n);
	vtype rdr=rot(rd2,rp1,rotinf.cc,rotinf.ss,rotinf.n);

	rotinftype rotinf180;
	vtype n=v(ez,cmr-rp1);n=(1/norm(n))*n;
	rotinf180.cc=-1;
	rotinf180.ss=0;
	rotinf180.n=n;

	vtype rdr180=rot(rdr,rp1,rotinf180.cc,rotinf180.ss,rotinf180.n);
	return norm(rdr180-rp2)>norm(rdr-rp2);
}


void clustertype::step2(const rotattempt& att,surroundtype& /*su*/,const vtype& rp)
{
	contacts=att.contacts;
	stable=att.stable;
	p1=att.p1;
	d1=att.d1;
	p2=att.p2;
	d2=att.d2;
	p3=att.p3;
	d3=att.d3;
	td=att.td;
	tp=att.tp;
	f=att.f;
	cm=rot(cm,rp,att.rotinf.cc,att.rotinf.ss,att.rotinf.n);
	for(int i=0;i<r.size();i++) r[i]=rot(r[i],rp,att.rotinf.cc,att.rotinf.ss,att.rotinf.n);
    if (orientation) *orientation=Quaternion::from_axis(att.rotinf.n,att.rotinf.cc,att.rotinf.ss)* *orientation;
}


void clustertype::find_steepest_descent(std::span<const vtype> rp,std::span<const vtype> rd,std::span<const int> p,std::span<const int> d,const vtype& cmc,rotattempt& att)
{
	gradrotattempt gatt;
	std::vector<gradrotattempt> gradevent;
	SDBD_REQUIRE(rp.size()==rd.size());
	SDBD_REQUIRE(p.size()==d.size());
	SDBD_REQUIRE(rp.size()==p.size());
	int contact_number=rp.size();
	//Add event for cluster falling, if allowed
	bool can_fall=true;
	for(int i=0;i<contact_number;i++) {can_fall=can_fall&&rd[i].z()<rp[i].z();}
	if(can_fall) 
	{
		gatt.contacts=0;
		gatt.stable=false;
		gatt.d1=-1;
		gatt.p1=-1;
		gatt.d2=-1;
		gatt.p2=-1;
		gatt.d3=-1;
		gatt.p3=-1;
		gatt.td=-1;
		gatt.tp=-1;
		gatt.f=0;
		gatt.cm=cmc;
		gatt.gz=-1;
		gradevent.push_back(gatt);
	}

	//Add trajectories for rolling on one contact, if they exist 
	for(int i=0;i<contact_number;i++)
	{
		bool can_rotate=true;
		for(int j=0;j<contact_number;j++) 
		if(i!=j) 
		{
			can_rotate=can_rotate&&engine_->g(rp[i],cmc,rd[j],rp[j])>engine_->epstang;
		}
		if(can_rotate)
		{
			gatt.contacts=1;
			gatt.stable=false;
			gatt.d1=d[i];
			gatt.p1=p[i];
			gatt.d2=-1;
			gatt.p2=-1;
			gatt.d3=-1;
			gatt.p3=-1;
			gatt.td=-1;
			gatt.tp=-1;
			gatt.f=0;
			gatt.cm=cmc;
			gatt.gz=s(ez,engine_->grad(rp[i],cmc));
			gradevent.push_back(gatt);
		}

	}

	//Add trajectories for rolling on two contacts, if they exist
	for(int i=1;i<contact_number;i++)
	for(int j=0;j<i;j++)
	{
		bool can_roll=true;
		for(int k=0;k<contact_number;k++) if(k!=i&&k!=j) {can_roll=can_roll&&engine_->g2(rp[i],rp[j],cmc,rd[k],rp[k])>engine_->epstang;}
		if(can_roll)
		{
			gatt.contacts=2;
			gatt.stable=false;
			gatt.d1=d[i];
			gatt.p1=p[i];
			gatt.d2=d[j];
			gatt.p2=p[j];
			gatt.d3=-1;
			gatt.p3=-1;
			gatt.td=-1;
			gatt.tp=-1;
			gatt.f=0;
			gatt.cm=cmc;
			gatt.gz=s(ez,engine_->grad2(rp[i],rp[j],cmc));
			gradevent.push_back(gatt);
		}
	}

	if(!gradevent.empty())//if there are trajectories find the steepest
	{
		gatt=gradevent[0];
		for(int ig=0;ig<gradevent.size();ig++)
		{
			if(gradevent[ig].gz<gatt.gz) gatt=gradevent[ig];
		}
		att.contacts=gatt.contacts;
		att.stable=gatt.stable;
		att.d1=gatt.d1;	
		att.p1=gatt.p1;
		att.d2=gatt.d2;
		att.p2=gatt.p2;
		att.d3=gatt.d3;
		att.p3=gatt.p3;
		att.td=gatt.td;
		att.tp=gatt.tp;
		att.f=gatt.f;
		att.cm=gatt.cm;
	}
	else//if there are no trajectories cluster is stable
	{
		att.contacts=rp.size();
		att.stable=true;
		att.d1=d[0];	
		att.p1=p[0];
		att.d2=d[1];
		att.p2=p[1];
		att.d3=d[2];
		att.p3=p[2];
		att.td=-1;
		att.tp=-1;
		att.f=0;
		att.cm=cmc;
	}
}

void clustertype::randrot()
{
	//random basis way

	vtype nxp;//first vector
	do
	{
		nxp.x()=engine_->ranfs();
		nxp.y()=engine_->ranfs();
		nxp.z()=engine_->ranfs();
	}
	while(norm(nxp)>1||norm(nxp)<0.5);
	nxp=(1/norm(nxp))*nxp;

	vtype nyp;//second vector
	do
	{
		nyp.x()=engine_->ranfs();
		nyp.y()=engine_->ranfs();
		nyp.z()=engine_->ranfs();
	}
	while(norm(nyp)>1||norm(nyp)<0.5||std::fabs(s(nyp,nxp))<0.1);

	nyp=nyp-nxp*s(nyp,nxp);//nyp is now normal to nxp
	nyp=(1/norm(nyp))*nyp;

	vtype nzp=v(nxp,nyp);//nzp normal to both nxp and nyp
    if (orientation) *orientation=Quaternion::from_basis({nxp,nyp,nzp})* *orientation;

	for(int i=0;i<r.size();i++) 
	{
		r[i]=cm+(r[i].x()-cm.x())*nxp+(r[i].y()-cm.y())*nyp+(r[i].z()-cm.z())*nzp;
	}
}
}
