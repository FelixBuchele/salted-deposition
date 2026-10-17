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

// Current or candidate contact state, with physical center of mass cm.
// Each d/p pair identifies a moving sphere and its fixed support sphere.
// Index -1 means absent; support -2 denotes the plane. td,tp remember a lost contact.
class clusterbase
{
	public:
		int contacts{};
		bool stable{};
		int d1{-1},p1{-1},d2{-1},p2{-1},d3{-1},p3{-1},td{-1},tp{-1};
        bool f{}; // Last contact change was a grazing loss requiring special recontact handling.
		vtype cm;
        // Compare future COM heights: during descent, the highest event is reached
        // first. max_element keeps the first candidate when heights tie.
		bool operator<(const clusterbase& c) const {return (*this).cm.z()<c.cm.z();}; 
        // Compare COM heights in descending order without changing event state.
		bool operator>(const clusterbase& c) const {return (*this).cm.z()>c.cm.z();}; 
};


// A possible end of vertical fall: future COM and contacts, without copying the spheres.
class fallattempt : public clusterbase
{

};


// A possible end of rotation: future COM, contacts and the rotation to reach them.
class rotattempt : public clusterbase
{
	public:
		rotinftype rotinf;
};


// A possible motion consistent with the contacts. gz is the vertical component
// of the unit COM direction; the most negative value gives the steepest descent.
class gradrotattempt : public rotattempt
{
	public:
		long double gz{};
};


// One moving particle: sphere centers r, radii R, contacts and physical COM.
// Sphere arrays use the same indices. Motion updates the stored COM rather than
// recalculating it from overlapping spheres. The engine must outlive this particle.
class clustertype : public clusterbase
{
	public:
		std::vector<vtype> r;
		std::vector<long double> R;

	public:
		// Calculate COM using radius-cubed weights, for the original built-in shapes.
		// Overlapping imported particles retain their supplied physical COM instead.
		vtype cmass();
		// Start a moving particle from validated geometry, retaining its physical COM,
		// particle id and exposed-sphere indices. Keep a reference to engine.
		clustertype(Engine& engine, const ParticleGeometry& geometry, std::size_t id);
		std::size_t particle_id;
        std::optional<Quaternion> orientation; // Engaged only for requested STL pose export.
		// Follow fall and rotation events until the particle stops in fixed environment su.
		// Return true if events stop making progress (numerical collapse); the caller
		// decides whether to accept that state or retry the drop.
		bool settle(surroundtype& su);
		// Translate body and COM together by a domain period when COM crosses a horizontal boundary.
		void periodic(long double lx,long double ly);
		// Translate all sphere centers and the physical COM by displacement v.
		// Do not recalculate COM from the spheres.
		void movecluster(const vtype& v);
		// Rotate constituents about the stored physical COM using the original random-basis
		// construction and GSL draw order.
		void randrot();
        // Shift geometry and COM horizontally by x,y in simulation length units.
        // The compatibility flag recalculates the built-in COM after shifting;
        // throw if it would overwrite an explicitly supplied physical COM.
        void shift_horizontal(long double x,long double y,bool legacy_recompute_center=false);
        // Return whether this sphere can contribute a contact surface.
        // Fully contained spheres are excluded; index must be within the body.
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
		// If the contacting sphere lies below its support rp, test a small phi rotation
		// before releasing contact. Append the event only if the height and side tests pass.
		void steeperdphi(const vtype& rp,std::vector<rotattempt>& atvec);
		// Search the swept rotation path about one support for new sphere contacts.
		// Overlapping search volumes cover the path; append events in search order.
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
		// Test a small phi rotation about supports rp1/rp2 when a steeper path would
		// use fewer contacts. Append it if the rotated state permits that path.
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
		// Choose the steepest allowed fall or rotation at COM cmc and write contact state att.
		// rp,rd are fixed/moving centers; p,d are their matching sphere indices.
		// Reject paths that close another contact; keep the first path when slopes tie.
		void find_steepest_descent(std::span<const vtype> rp,std::span<const vtype> rd,std::span<const int> p,std::span<const int> d,const vtype& cmc,rotattempt& att);

};
}
