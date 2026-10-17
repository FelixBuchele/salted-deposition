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
// One sphere of a rigid particle, in body coordinates and the particle's declared units.
struct Sphere {
    Vec3 center{};
    long double radius{};
};

// Axis-aligned minimum/maximum coordinates of the sphere union, including radii.
struct Bounds {
    Vec3 minimum;
    Vec3 maximum;
};

// Physical properties of the target geometry, assuming uniform density.
// Inertia is about the physical center of mass (COM), in the current axes and
// at unit density (length^5). Overlapping collision spheres do not define mass.
struct PhysicalProperties {
    long double volume{};
    long double density{1};
    std::array<long double,9> inertia_per_density{}; // Row-major symmetric tensor.
    std::string source, checksum, backend_version, parameters, source_units;
    long double source_scale{1}, voxel_size{};
    // Return density times target volume, independent of collision-sphere overlaps.
    long double mass() const { return density*volume; }
};

// Location and orientation of the original STL coordinate frame within the particle.
// Its uniform scale is PhysicalProperties::source_scale; the mesh itself is not copied.
struct StlSource {
    std::filesystem::path path;
    Vec3 origin;
    Quaternion orientation;
};

// Reusable particle geometry, separate from the deposition state.
// For disjoint spheres, COM may be calculated at uniform density. For overlapping
// spheres, the caller must supply the physical COM: sphere weighting would count
// shared volume more than once. Validation cannot establish its physical accuracy.
class ParticleGeometry {
public:
    // Bound initial pairwise validation to fewer than 50 million sphere pairs.
    static constexpr std::size_t maximum_spheres = 10'000;
    // Store and validate spheres, COM and optional physical properties.
    // Overlaps require an explicit COM. The units label does not rescale coordinates.
    // Throw for unknown units, invalid geometry or values outside the supported range.
    explicit ParticleGeometry(std::vector<Sphere> spheres,
                              std::optional<Vec3> physical_center = std::nullopt,
                              std::optional<PhysicalProperties> properties = std::nullopt,
                              std::optional<StlSource> source = std::nullopt,
                              std::string units = "legacy");
    // Return the declared coordinate unit; legacy means caller-assigned original units.
    const std::string& units() const noexcept { return units_; }
    // Return the original STL reference and frame, if available for pose export.
    const std::optional<StlSource>& stl_source() const noexcept { return source_; }
    // View the spheres in declared units; valid until this geometry is replaced or destroyed.
    std::span<const Sphere> spheres() const noexcept { return spheres_; }
    // Return the stored physical COM in the same coordinates as constituent centers.
    Vec3 center_of_mass() const noexcept { return center_; }
    // Return optional target volume, density, inertia about COM and source information.
    // The reference remains valid while this object is alive.
    const std::optional<PhysicalProperties>& physical_properties() const noexcept { return properties_; }
    // Return axis-aligned bounds including radii, in declared length units.
    Bounds bounds() const noexcept { return bounds_; }
    // Return the bounding-box dimensions including radii, in declared length units.
    Vec3 dimensions() const noexcept { return bounds_.maximum - bounds_.minimum; }
    // Return whether the stored COM was explicitly supplied or carried through a transformation.
    bool has_explicit_center() const noexcept { return explicit_center_; }
    // Return whether any constituent pair strictly overlaps; internal overlaps belong to
    // the same rigid body.
    bool has_overlaps() const noexcept { return overlaps_; }
    // Return whether the collision union is a single sphere centered exactly on the physical COM.
    bool is_centered_sphere() const noexcept { return centered_sphere_; }
    // Return the largest sphere radius, in declared length units.
    long double maximum_radius() const noexcept { return maximum_radius_; }
    // Return the radius of a COM-centered bound enclosing all spheres, in declared units.
    long double bounding_radius() const noexcept { return bounding_radius_; }
    // Return indices of spheres not fully contained in another sphere.
    // Contained spheres remain stored and exported, but cannot add a contact surface.
    std::span<const std::size_t> collision_spheres() const noexcept { return collision_spheres_; }

    // Translate spheres, COM and the STL frame by displacement; return a new geometry.
    // Inertia about COM is unchanged. Throw if the resulting coordinates are invalid.
    ParticleGeometry translated(Vec3 displacement) const;
    // Scale all lengths about the origin by a positive factor; return a new geometry.
    // Volume scales by factor^3 and unit-density inertia by factor^5. Throw on invalid values.
    ParticleGeometry scaled(long double factor) const;
    // Rotate spheres, COM and the STL frame about axis through pivot; return a new geometry.
    // The angle is in radians. Rotate inertia as R I R^T; throw on invalid inputs or results.
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
// controls discretization; random supplies uniform numbers in [0,1). length_scale
// scales radii, spacing and jitter, leaving the dropping origin unchanged.
ParticleGeometry generate_builtin(BuiltinShape shape, int resolution,
                                  const UniformRandom& random, long double initial_height = 0,
                                  long double length_scale = 1);

// Return a symmetric stick of equal, overlapping spheres along x. Length includes
// end caps; radius and length must use the same unit. Symmetry gives COM=0 at uniform density.
// Reject nonfinite/nonpositive dimensions, fewer than two spheres, excessive count,
// or spacing outside (0,2*radius). No random numbers are consumed.
ParticleGeometry generate_stick(long double length, long double radius, std::size_t sphere_count);

// Read version-1/2/3 particle data using a decimal point; apply its scale and validate.
// Resolve version-3 STL paths relative to base (current directory if omitted).
// Return geometry with its physical COM, or throw on invalid input.
ParticleGeometry read_particle(std::istream& input,const std::filesystem::path& base = {});
// Write geometry, COM, optional properties and the STL frame with enough digits
// to recover long-double values. Write STL paths relative to base when supplied.
// Throw on output failure.
void write_particle(std::ostream& output, const ParticleGeometry& particle,const std::filesystem::path& base = {});
// Load and validate a particle file; return geometry or throw with path-qualified parsing
// diagnostics.
ParticleGeometry load_particle(const std::filesystem::path& path);
// Save the particle in the format above, replacing an existing file.
// Throw if opening, writing or closing the file fails.
void save_particle(const std::filesystem::path& path, const ParticleGeometry& particle);
}
