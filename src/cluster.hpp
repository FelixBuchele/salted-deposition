// Adapted from the original SDBD source by Nikola Topic (2015).
// See Topic & Pöschel, JCP 308 (2016), DOI: 10.1016/j.jcp.2015.12.052.
#pragma once
#include "math/rotation.hpp"
#include <span>
#include "surround.hpp"
#include "sdbd/particle.hpp"
namespace sdbd::detail {
class surroundtype;
class Engine;

// Contact-event state and physical COM; d indices name moving constituents, p indices fixed
// spheres; -1 denotes absent contacts.
class clusterbase
{
	public:
		int contacts{};
		bool stable{};
		int d1{-1},p1{-1},d2{-1},p2{-1},d3{-1},p3{-1},td{-1},tp{-1};
        bool f{}; // Last contact change was a grazing loss requiring special recontact handling.
		vtype cm;
        // Compare future COM heights; max_element retains the first equal-height event.
		bool operator<(const clusterbase& c) const {return (*this).cm.z()<c.cm.z();}; 
        // Compare COM heights in descending order without changing event state.
		bool operator>(const clusterbase& c) const {return (*this).cm.z()>c.cm.z();}; 
};


// Candidate falling event containing contact state and future physical COM, without constituent
// copies.
class fallattempt : public clusterbase
{

};


// Candidate rotational event with contact state, physical COM and common axis/cosine/sine data.
class rotattempt : public clusterbase
{
	public:
		rotinftype rotinf;
};


// Candidate contact-constrained trajectory; gz is the vertical component of its unit COM
// descent direction.
class gradrotattempt : public rotattempt
{
	public:
		long double gz{};
};


// Mutable rigid-particle solver state borrowing one Engine; coordinates, radii and collision
// flags share indexing and physical COM is authoritative.
class clustertype : public clusterbase
{
	public:
		std::vector<vtype> r;
		std::vector<long double> R;

	public:
		// Return the sphere-volume-weighted COM of disjoint spheres r,R in input order; do not use
		// for overlapping-particle physics.
		vtype cmass();
		// Copy validated constituent geometry and physical COM into mutable body state; borrow
		// engine and preserve particle id/collision membership.
		clustertype(Engine& engine, const ParticleGeometry& geometry, std::size_t id);
		std::size_t particle_id;
        std::optional<Quaternion> orientation; // Engaged only for requested STL pose export.
		// Advance contact events against su until fixed; return true for a collapsed attempt that
		// may require redrop, preserving physical COM.
		bool settle(surroundtype& su);
		// Translate body and COM together by a domain period when COM crosses a horizontal boundary.
		void periodic(long double lx,long double ly);
		// Translate every constituent and physical COM by v in legacy length units without
		// recomputing mass distribution.
		void movecluster(const vtype& v);
		// Rotate constituents about the stored physical COM using the original random-basis
		// construction and GSL draw order.
		void randrot();
        // Translate x/y geometry and COM; legacy_recompute_center is allowed only for
        // inferred disjoint-sphere COM, otherwise throw.
        void shift_horizontal(long double x,long double y,bool legacy_recompute_center=false);
        // Return whether constituent index contributes an exposed union boundary; index
        // must lie within the body.
        bool collision_active(std::size_t index) const { return collision_active_[index]; }

	private:
        bool explicit_center_{};
        Engine* engine_; // non-owning; Engine::execute owns the particle vector
        std::vector<bool> collision_active_;
		// Choose the highest future falling event from fixed spheres and plane; update body or
		// throw when no event exists.
		void fall(surroundtype& su);
		// Append future sphere-contact events for vertical fall using periodic grid candidates and
		// exposed constituents.
		void bfall(surroundtype& su,std::vector<fallattempt>& attempts);
		// Append plane-contact events for exposed constituents; the original model fixes the body
		// on first contact with z=0.
		void gfall(std::vector<fallattempt>& attempts);
		// Apply falling event att and its contact state; su supplies the fixed support for the
		// legacy contact-position correction.
		void step0(fallattempt att,surroundtype& su);
		// Move body/COM downward by dlt; correct only constituent d1 against support p1 for legacy
		// cancellation error, allowing small rigidity drift.
		void translate(long double dlt,int d1,int p1,surroundtype& su);

		// Choose/apply the next one-contact event in su; temporarily perturb a directly-above
		// support with two horizontal random draws.
		void rotate1(surroundtype& su);
		// Append permitted phi-rotation events for a frustrated one-support trajectory about rp;
		// retain original height/side tests.
		void steeperdphi(const vtype& rp,std::vector<rotattempt>& atvec);
		// Append collision events found along one-support rotation using overlapping grid-search
		// volumes; preserve candidate and tie order.
		void brot1(surroundtype& s,std::vector<rotattempt>& attempts);
		// Append loss-of-contact event when the contacting constituent reaches support rp height.
		void hrot(const vtype& rp,std::vector<rotattempt>& attempts);
		// Append rotations about rp contacting the plane; plane contact fixes the body.
		void grot1(const vtype& rp,std::vector<rotattempt>& attempts);
		// Append stable hanging event with physical COM below support rp when geometrically permitted.
		void bottomrot1(const vtype& rp,std::vector<rotattempt>&  attempts);
		// Apply event contact state and rotate every constituent and physical COM about rp with the
		// same rotation.
		void step1(const rotattempt& att,const vtype& rp);

		// Choose/apply the next two-contact event in su, including grazing loss, frustration,
		// collision, hanging and plane contact.
		void rotate2(surroundtype& su);
		// Append phi-rotation events that resolve a frustrated trajectory about supports rp1/rp2 by
		// selecting a steeper valid contact set.
		void steeperdphi2(const vtype& rp1,const vtype& rp2,std::vector<rotattempt>& atvec);
		// Append sphere-contact events during two-support rotation using periodic grid candidates
		// and unchanged event order.
		void brot2(surroundtype& su,std::vector<rotattempt>&  attempts);
		// Append stable hanging event with physical COM below the rp1/rp2 support axis.
		void bottomrot2(const vtype& rp1,const vtype& rp2,std::vector<rotattempt>& atvec);
		// Append plane-contact rotations about the rp1/rp2 axis; plane contact fixes the body.
		void grot2(const vtype rp1,const vtype rp2,std::vector<rotattempt>&  attempts);
		// Append simultaneous loss-of-two-contacts events when a constituent reaches its support height.
		void hrot2(const vtype& rp1,const vtype& rp2,std::vector<rotattempt>& atvec);
		// Append permitted grazing contact-loss events about rp1/rp2 using paper Eqs. 25 and 31; su
		// is unused.
		void disconnect(const vtype& rp1,const vtype& rp2, std::vector<rotattempt>&  attempts,surroundtype& su);
		// Test the sign condition in paper Eq. 31 by comparing rd2-support distance before/after a
		// half turn; return whether it increases.
		bool positivecircle(const vtype& rp1,const vtype& rp2,const vtype& rd2,const rotinftype& rotinf);
		// Apply event state and identical rotation to constituents/COM about rp; su is unused.
		void step2(const rotattempt& att,surroundtype& su,const vtype& rp);
		// Select the valid falling/one-/two-support path with lowest vertical gradient; read
		// contact spans and write att, retaining tie order.
		void find_steepest_descent(std::span<const vtype> rp,std::span<const vtype> rd,std::span<const int> p,std::span<const int> d,const vtype& cmc,rotattempt& att);

};
}
