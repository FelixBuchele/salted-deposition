#pragma once
#include "sdbd/vector.hpp"
#include "sdbd/orientation.hpp"
#include <array>
#include <string>
#include <filesystem>
#include <functional>
#include <iosfwd>
#include <optional>
#include <span>
#include <vector>

namespace sdbd {
// One constituent sphere in body coordinates; center and positive radius use legacy length units.
struct Sphere {
    Vec3 center{};
    long double radius{};
};

// Axis-aligned minimum/maximum coordinates of the sphere union, including radii.
struct Bounds {
    Vec3 minimum;
    Vec3 maximum;
};

// Uniform density. Inertia is about the physical COM in the current axes,
// at unit density (length^5); mass = density * target volume.
struct PhysicalProperties {
    long double volume{};
    long double density{1};
    std::array<long double,9> inertia_per_density{}; // Row-major symmetric tensor.
    std::string source, checksum, backend_version, parameters, source_units;
    long double source_scale{1}, voxel_size{};
    // Return density times target volume, independent of collision-sphere overlaps.
    long double mass() const { return density*volume; }
};

// Original STL reference and its current body-frame origin/active rotation.
// Uniform scale is PhysicalProperties::source_scale; no mesh data is stored here.
struct StlSource {
    std::filesystem::path path;
    Vec3 origin;
    Quaternion orientation;
};

// Geometry only: no contact state, world position, RNG ownership or solver dependency.
// Automatic COM assumes uniform density and disjoint spheres. Overlaps require
// a caller-supplied physical COM; its physical correctness cannot be inferred.
class ParticleGeometry {
public:
    // Bound initial pairwise validation to fewer than 50 million sphere pairs.
    static constexpr std::size_t maximum_spheres = 10'000;
    // Validate finite spheres and optional physical COM/properties; overlaps require
    // explicit COM. Units label coordinates without converting them; reject unknown
    // units, invalid geometry or numerical range.
    explicit ParticleGeometry(std::vector<Sphere> spheres,
                              std::optional<Vec3> physical_center = std::nullopt,
                              std::optional<PhysicalProperties> properties = std::nullopt,
                              std::optional<StlSource> source = std::nullopt,
                              std::string units = "legacy");
    // Return the declared coordinate unit; legacy means caller-assigned original units.
    const std::string& units() const noexcept { return units_; }
    // Return original STL/frame provenance when retained; absence disables STL-pose export.
    const std::optional<StlSource>& stl_source() const noexcept { return source_; }
    // Return a read-only constituent view in legacy length units; lifetime follows this object.
    std::span<const Sphere> spheres() const noexcept { return spheres_; }
    // Return the stored physical COM in the same coordinates as constituent centers.
    Vec3 center_of_mass() const noexcept { return center_; }
    // Return optional target-derived volume, density, COM-frame inertia and provenance;
    // lifetime follows this geometry.
    const std::optional<PhysicalProperties>& physical_properties() const noexcept { return properties_; }
    // Return axis-aligned bounds including sphere radii, in legacy length units.
    Bounds bounds() const noexcept { return bounds_; }
    // Return positive axis-aligned extent including radii, in legacy length units.
    Vec3 dimensions() const noexcept { return bounds_.maximum - bounds_.minimum; }
    // Return whether the stored COM was explicitly supplied or carried through a transformation.
    bool has_explicit_center() const noexcept { return explicit_center_; }
    // Return whether any constituent pair strictly overlaps; internal overlaps belong to
    // the same rigid body.
    bool has_overlaps() const noexcept { return overlaps_; }
    // Return whether the collision union is a single sphere centered exactly on the physical COM.
    bool is_centered_sphere() const noexcept { return centered_sphere_; }
    // Return largest constituent radius in legacy length units.
    long double maximum_radius() const noexcept { return maximum_radius_; }
    // Return the maximum COM-to-center distance plus radius, in legacy length units.
    long double bounding_radius() const noexcept { return bounding_radius_; }
    // Return exposed constituent indices; contained spheres remain stored for I/O.
    std::span<const std::size_t> collision_spheres() const noexcept { return collision_spheres_; }

    // Return geometry and physical COM translated by displacement; inertia about COM is
    // unchanged. Reject nonfinite/range-invalid results.
    ParticleGeometry translated(Vec3 displacement) const;
    // Return geometry scaled about the origin by positive factor; volume scales cubically
    // and unit-density inertia by factor^5. Reject invalid range.
    ParticleGeometry scaled(long double factor) const;
    // Return geometry/COM rotated by radians about nonzero axis through pivot; rotate
    // COM-frame inertia by R I R^T. Reject invalid inputs/range.
    ParticleGeometry rotated(Vec3 axis, long double radians, Vec3 pivot = Vec3{}) const;

private:
    std::string units_;
    // Recompute finite coordinate bounds and radii in O(number of spheres).
    void update_bounds();
    // Validate COM bounds and update the maximum COM-to-surface distance.
    void update_center_bounds();
    // Reject invalid target volume, density, tensor or provenance data.
    void validate_properties() const;
    std::optional<PhysicalProperties> properties_;
    std::optional<StlSource> source_;
    std::vector<Sphere> spheres_;
    Vec3 center_;
    Bounds bounds_;
    bool explicit_center_{};
    bool overlaps_{};
    bool centered_sphere_{};
    long double maximum_radius_{};
    long double bounding_radius_{};
    std::vector<std::size_t> collision_spheres_;
};

enum class BuiltinShape { sphere, spiky, cube };
// Callback returns one original uniform variate in [0,1). No RNG substitution.
using UniformRandom = std::function<long double()>;
// Generate original sphere, six-arm or surface-lattice cube at initial_height; resolution
// controls discretization and random supplies baseline variates. length_scale
// scales radii, spacing and jitter, leaving the dropping origin unchanged.
ParticleGeometry generate_builtin(BuiltinShape shape, int resolution,
                                  const UniformRandom& random, long double initial_height = 0,
                                  long double length_scale = 1);

// Return a symmetric stick of equal, overlapping spheres along x. Length includes
// end caps; radius and length use legacy units. Uniform union density gives COM=0.
// Reject nonfinite/nonpositive dimensions, fewer than two spheres, excessive count,
// or spacing outside (0,2*radius). No random numbers are consumed.
ParticleGeometry generate_stick(long double length, long double radius, std::size_t sphere_count);

// Parse version-1/2/3 multisphere data in classic locale; return validated scaled
// geometry/COM or throw. Resolve version-3 STL paths against base (cwd if omitted).
ParticleGeometry read_particle(std::istream& input,const std::filesystem::path& base = {});
// Write canonical max_digits10 geometry, COM, optional properties and STL frame;
// express STL paths relative to nonempty base. Throw on output failure.
void write_particle(std::ostream& output, const ParticleGeometry& particle,const std::filesystem::path& base = {});
// Load and validate a particle file; return geometry or throw with path-qualified parsing
// diagnostics.
ParticleGeometry load_particle(const std::filesystem::path& path);
// Write canonical particle geometry to path, replacing an existing file; throw on
// open/write/close failure.
void save_particle(const std::filesystem::path& path, const ParticleGeometry& particle);
}
