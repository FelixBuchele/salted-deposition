# Particle files, properties and generation

All representations use Cartesian body coordinates with finite centers and
strictly positive radii. Coordinates declare mm, cm, m or unassigned legacy units. The simulation
uses the same declared unit; legacy files are interpreted in the run unit, with
no implicit rescaling. Particle bounds, COM, numerical range and domain compatibility are
validated. Internal overlaps are permitted. Automatic COM is available only for
disjoint manual spheres; it is the sphere-volume-weighted COM at uniform density.
Overlapping manual geometry requires `center explicit x y z`.

## File format

Files are ordered whitespace-delimited records in classic locale. The complete
minimal example below illustrates an overlapping manual definition (unassigned legacy units):

```text
SDBD_MULTISPHERE 1
units legacy
scale 1
center explicit 0 0 0
spheres 2
-15 0 0 30
15 0 0 30
```

`scale` is a positive factor applied to coordinates, radii and COM while loading.
Canonical saves write `scale 1` and represented values at `max_digits10` precision.
The `units` record can be `legacy`, `mm`, `cm` or `m`; this label is preserved
through transformations and I/O and does not itself convert coordinates.
`center auto` selects the disjoint-sphere approximation. Exactly the declared
number of x/y/z/r rows must follow `spheres`; trailing data is rejected. Unknown
versions, malformed/nonfinite numbers, invalid radii and excess counts fail
explicitly. Comments are not part of this particle-file grammar.

Version 2 inserts these records before `spheres`, in this order:

| Record | Meaning |
|---|---|
| `target_volume`, `density` | Positive target volume in L³ and uniform mass/L³ |
| `inertia_per_density` | Nine row-major symmetric COM-frame tensor entries in L⁵ |
| `source`, `checksum`, `backend`, `parameters`, `source_units` | Quoted provenance strings |
| `source_scale`, `voxel_size` | Original-source-to-solver factor and target voxel spacing in L |

Version 3 additionally inserts `source_stl` (quoted original mesh path),
`source_origin` (three coordinates) and `source_orientation_wxyz` (unit quaternion)
after the version-2 records. Paths are resolved relative to the particle file.
This retains the original STL frame through transforms and I/O. Physical COM is
stored separately in `center`; it is not inferred from the source origin.
Versions 1/2 remain loadable but lack authoritative original-mesh frame data for
STL pose export. See [export conventions](exports.md).

## Generators

Built-ins preserve Topic's sphere, six-arm spiky representation and perturbed
cube surface lattice. The cube at resolution 3 has 26 spheres, not a filled
27-sphere volume; its legacy COM approximation is retained for comparisons.
Runtime generation consumes the original GSL stream in reference order.

`particle = stick` generates overlapping equal spheres along x with an explicit COM at
the origin. Length includes end caps; spacing is (length−2r)/(n−1), with
0 < spacing < 2r. Uniform union density and symmetry justify this COM.

The optional GEMSS adapter accepts binary STL with explicit source units/scale.
It checks finite triangles, nondegeneracy, closure, orientation, connectivity and
resource bounds where supported. Target-grid COM, volume and inertia are mandatory;
failure is an error, never a fallback to sphere weighting. `maximum_spheres`
controls collision resolution; `voxel_resolution` controls the target grid.
Physical properties at fixed voxel resolution do not change with sphere count.
Store both achieved and requested reconstruction accuracy and inspect the result.
See [configuration reference](configuration.md) and [rock example](../examples/rock/README.md) for commands.
