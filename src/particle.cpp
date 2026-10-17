#include "sdbd/particle.hpp"
#include "geometry_validation.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace sdbd {
namespace {
using geometry_detail::ValidationEnvironment;
using geometry_detail::arithmetic_range;
// Return whether all Cartesian components of p are finite.
bool finite(Vec3 p) { return std::isfinite(p.x()) && std::isfinite(p.y()) && std::isfinite(p.z()); }
// Throw invalid_argument with message when a geometry invariant is false.
void require(bool valid, const char* message) {
    if (!valid) throw std::invalid_argument(message);
}
// Return robust Euclidean distance between a,b in their shared length units.
long double distance(Vec3 a, Vec3 b) {
    const Vec3 delta = a-b;
    return std::hypot(delta.x(), delta.y(), delta.z());
}
}

ParticleGeometry::ParticleGeometry(std::vector<Sphere> spheres, std::optional<Vec3> physical_center,
                                   std::optional<PhysicalProperties> properties,std::optional<StlSource> source,std::string units)
    : units_(std::move(units)), properties_(std::move(properties)), source_(std::move(source)), spheres_(std::move(spheres)), explicit_center_(physical_center.has_value()) {
    ValidationEnvironment environment;
    require(units_=="legacy" || units_=="mm" || units_=="cm" || units_=="m", "Unknown particle coordinate unit");
    require(!spheres_.empty() && spheres_.size() <= maximum_spheres, "Sphere count must be between 1 and 10000");
    update_bounds();
    std::vector<bool> exposed(spheres_.size(),true);
    for (std::size_t i=0; i<spheres_.size(); ++i) {
        for (std::size_t j=0; j<i; ++j) {
            const auto d = distance(spheres_[i].center, spheres_[j].center);
            const auto radii = spheres_[i].radius+spheres_[j].radius;
            require(std::isfinite(d) && std::isfinite(radii), "Particle distance exceeds numerical range");
            if (d<radii) overlaps_=true;
            // Exact containment removes no union boundary. Keep the first of
            // coincident equal spheres; retain all spheres in geometry/output.
            if (d+spheres_[i].radius<=spheres_[j].radius) exposed[i]=false;
            else if (d+spheres_[j].radius<=spheres_[i].radius) exposed[j]=false;
        }
    }
    for (std::size_t i=0;i<exposed.size();++i) if (exposed[i]) collision_spheres_.push_back(i);
    if (physical_center) {
        center_=*physical_center;
    } else {
        require(!overlaps_, "Overlapping spheres require an explicit physical center of mass");
        // Preserve the original summation order and expression precision.
        long double mass=0;
        Vec3 sum;
        for (const auto& sphere : spheres_) {
            const long double r3=sphere.radius*sphere.radius*sphere.radius;
            require(std::isfinite(r3) && r3>0, "Sphere volume exceeds numerical range");
            sum=sum+r3*(sphere.center-spheres_[0].center);
            mass=mass+r3;
            require(finite(sum) && std::isfinite(mass), "Center of mass exceeds numerical range");
        }
        center_=(1/mass)*sum;
        center_=center_+spheres_[0].center;
    }
    update_center_bounds();
    // Exact union test, not an approximate shape simplification. Keep every
    // constituent sphere; only select the existing centered-sphere collapse rule.
    for (const auto& outer : spheres_) {
        if (outer.center.x()!=center_.x() || outer.center.y()!=center_.y() || outer.center.z()!=center_.z()) continue;
        const bool contains_all=std::all_of(spheres_.begin(),spheres_.end(),[&](const Sphere& inner) {
            return distance(outer.center,inner.center)+inner.radius<=outer.radius;
        });
        if (contains_all) { centered_sphere_=true; break; }
    }
    validate_properties();
    arithmetic_range();
}

void ParticleGeometry::update_bounds() {
    maximum_radius_=0;
    const long double infinity = std::numeric_limits<long double>::infinity();
    bounds_ = {Vec3(infinity,infinity,infinity), Vec3(-infinity,-infinity,-infinity)};
    for (const auto& sphere : spheres_) {
        require(finite(sphere.center) && std::isfinite(sphere.radius) && sphere.radius > 0,
                "Sphere coordinates must be finite and radius positive and finite");
        const auto p = sphere.center;
        const auto r = sphere.radius;
        bounds_.minimum = Vec3(std::min(bounds_.minimum.x(),p.x()-r), std::min(bounds_.minimum.y(),p.y()-r), std::min(bounds_.minimum.z(),p.z()-r));
        bounds_.maximum = Vec3(std::max(bounds_.maximum.x(),p.x()+r), std::max(bounds_.maximum.y(),p.y()+r), std::max(bounds_.maximum.z(),p.z()+r));
        maximum_radius_ = std::max(maximum_radius_,r);
    }
    const auto size = dimensions();
    require(finite(bounds_.minimum) && finite(bounds_.maximum) && finite(size) &&
            size.x()>0 && size.y()>0 && size.z()>0, "Invalid particle dimensions or bounds");
}

void ParticleGeometry::update_center_bounds() {
    bounding_radius_=0;
    require(finite(center_) && center_.x()>=bounds_.minimum.x() && center_.x()<=bounds_.maximum.x() &&
            center_.y()>=bounds_.minimum.y() && center_.y()<=bounds_.maximum.y() &&
            center_.z()>=bounds_.minimum.z() && center_.z()<=bounds_.maximum.z(),
            "Physical center of mass is nonfinite or outside particle bounds");
    for (const auto& sphere : spheres_)
        bounding_radius_=std::max(bounding_radius_,distance(sphere.center,center_)+sphere.radius);
    require(std::isfinite(bounding_radius_) && bounding_radius_>0, "Invalid particle bounding radius");
}

void ParticleGeometry::validate_properties() const {
    if (source_) {
        require(properties_.has_value() && !source_->path.empty() && finite(source_->origin),"Invalid STL source frame");
        const auto& q=source_->orientation;
        const auto length=std::hypot(std::hypot(q.w,q.x),std::hypot(q.y,q.z));
        require(std::isfinite(length) && std::abs(length-1)<1e-12L,"STL source orientation must be a unit quaternion");
    }
    if (properties_) {
        require(explicit_center_, "Physical properties require an explicit physical COM");
        const auto& p=*properties_;
        require(std::isfinite(p.volume) && p.volume>0 && std::isfinite(p.density) && p.density>0 &&
                std::isfinite(p.mass()) && p.mass()>0, "Invalid physical volume or density");
        require(std::isfinite(p.source_scale) && p.source_scale>0 && std::isfinite(p.voxel_size) && p.voxel_size>0,
                "Invalid physical-property scale or voxel size");
        const auto& a=p.inertia_per_density;
        for (auto value:a) require(std::isfinite(value), "Nonfinite inertia tensor");
        require(a[1]==a[3] && a[2]==a[6] && a[5]==a[7], "Inertia tensor must be symmetric");
        // Sylvester criterion: reject nonphysical/degenerate volumetric tensors.
        require(a[0]>0 && a[0]*a[4]-a[1]*a[1]>0 &&
                a[0]*(a[4]*a[8]-a[5]*a[7])-a[1]*(a[3]*a[8]-a[5]*a[6])+a[2]*(a[3]*a[7]-a[4]*a[6])>0,
                "Inertia tensor must be positive definite");
        require(!p.source.empty() && !p.checksum.empty() && !p.backend_version.empty() &&
                !p.parameters.empty() && !p.source_units.empty(), "Incomplete physical-property provenance");
    }
}

ParticleGeometry ParticleGeometry::translated(Vec3 displacement) const {
    ValidationEnvironment environment;
    require(finite(displacement), "Translation must be finite");
    auto result=*this;
    auto& spheres=result.spheres_;
    for (auto& sphere : spheres) sphere.center+=displacement;
    if (result.source_) result.source_->origin+=displacement;
    const auto center=center_+displacement;
    arithmetic_range();
    // A rigid transform moves the physical COM, even when initially inferred.
    result.center_=center;result.explicit_center_=true;
    result.update_bounds();result.update_center_bounds();
    arithmetic_range();
    return result;
}

ParticleGeometry ParticleGeometry::scaled(long double factor) const {
    ValidationEnvironment environment;
    require(std::isfinite(factor) && factor>0, "Scale must be positive and finite");
    auto result=*this;
    auto& spheres=result.spheres_;
    for (auto& sphere : spheres) { sphere.center*=factor; sphere.radius*=factor; }
    if (result.source_) result.source_->origin*=factor;
    const auto center=factor*center_;
    auto properties=properties_;
    if (properties) {
        properties->volume*=factor*factor*factor;
        for (auto& value:properties->inertia_per_density) value*=factor*factor*factor*factor*factor;
        properties->source_scale*=factor;
        properties->voxel_size*=factor;
    }
    arithmetic_range();
    result.center_=center;result.explicit_center_=true;result.properties_=std::move(properties);
    result.update_bounds();result.update_center_bounds();result.validate_properties();
    arithmetic_range();
    return result;
}

ParticleGeometry ParticleGeometry::rotated(Vec3 axis, long double radians, Vec3 pivot) const {
    ValidationEnvironment environment;
    require(finite(axis) && finite(pivot) && std::isfinite(radians), "Rotation data must be finite");
    const long double length=std::hypot(axis.x(),axis.y(),axis.z());
    require(std::isfinite(length) && length>0, "Rotation axis must be nonzero");
    axis=(1/length)*axis;
    const auto cc=std::cos(radians), ss=std::sin(radians);
    auto result=*this;
    auto& spheres=result.spheres_;
    for (auto& sphere : spheres) sphere.center=rot(sphere.center,pivot,cc,ss,axis);
    if (result.source_) {
        result.source_->origin=rot(result.source_->origin,pivot,cc,ss,axis);
        result.source_->orientation=Quaternion::from_axis(axis,cc,ss)*result.source_->orientation;
    }
    const auto center=rot(center_,pivot,cc,ss,axis);
    auto properties=properties_;
    if (properties) {
        const Vec3 columns[]={rot(ex,Vec3{},cc,ss,axis),rot(ey,Vec3{},cc,ss,axis),rot(ez,Vec3{},cc,ss,axis)};
        long double r[3][3];
        for (int j=0;j<3;++j) { r[0][j]=columns[j].x(); r[1][j]=columns[j].y(); r[2][j]=columns[j].z(); }
        const auto old=properties->inertia_per_density;
        for (int i=0;i<3;++i) for (int j=i;j<3;++j) {
            long double value=0;
            for (int k=0;k<3;++k) for (int l=0;l<3;++l) value+=r[i][k]*old[3*k+l]*r[j][l];
            properties->inertia_per_density[3*i+j]=properties->inertia_per_density[3*j+i]=value;
        }
    }
    arithmetic_range();
    result.center_=center;result.explicit_center_=true;result.properties_=std::move(properties);
    result.update_bounds();result.update_center_bounds();result.validate_properties();
    arithmetic_range();
    return result;
}
}
