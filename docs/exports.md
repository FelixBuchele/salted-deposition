# Final deposition exports

Exports use only C++20/standard-library streams. The existing
GSL dependency is unchanged; neither VTK nor ParaView is linked into SALTED.
Every normal run writes sphere VTK, and prepared STL runs additionally write
poses when retained source frames are available. Set `stl_poses = false` for
collision-only output. The saved run record retains that choice. For example:

```sh
build/deposit run examples/rock/config.conf
```

Prepare the STL once using the [rock generation configuration](../examples/rock/generation.conf).
Generation writes a reusable particle file, not a deposition dataset. Default run
outputs are only `run.conf`, `deposition.vtk`, and applicable `stl_poses.vtk`.
The compatibility numerical files are optional through `legacy_output = true`.

| Optional file | Contents |
|---|---|
| `positions.txt` | Wrapped constituent centers and radii as x/y/z/r rows |
| `height.txt` | Twice the packing COM height; not the maximum surface height |
| `details.txt` | Particle count and reported interparticle intersections |
| `particles.txt` | Particle identifiers, physical COM coordinates and redrop counts |

`run.conf` records the actual maximum surface height as `surface_height`.

## Multisphere collision geometry

`deposition.vtk` is ASCII legacy VTK PolyData: one point/vertex cell per stored
constituent, with point-associated arrays `radius` (one component), `particle_id`
and `sphere_id` (zero-based integers). IDs are preserved even when spheres overlap
within one body. Contained constituents are retained. Coordinates are the solver's
final stored world coordinates, including its horizontal periodic wrapping.

All arrays use a `FIELD PointAttributes` block following `POINT_DATA`. This avoids
legacy readers silently loading only the first scalar array. Decimal coordinates,
radii and poses use `long double::max_digits10`; VTK's portable numeric type is
`double`. Text preserves solver digits, but ParaView reads binary64 and its GPU
instance transforms can use binary32. A rendered view is not a long-double
numerical verification.

In ParaView, open `deposition.vtk`, apply the reader and select **3D Glyphs** as
the representation (not the Glyph filter). In advanced display properties choose
Sphere, source radius **0.5**, enable scaling, select `radius`, use Magnitude
scaling and global Scale Factor **2**. Disable orientation, masking and glyph
culling/LOD for a full-resolution review. Color by point `particle_id`. Reset the
camera and zoom out if glyph extents exceed the initially fitted point bounds.
No triangulated sphere surfaces are stored in the file.

## Original STL surface instances

`stl_poses.vtk` contains one point/vertex per accepted rigid particle. Its point
arrays are `particle_id`, `orientation` (four-component unit quaternion, **WXYZ**)
and `scale` (one positive uniform factor). Quaternion rotations are right-handed,
active, and applied to original STL coordinates. Updates compose on the left:
the latest world rotation acts after the previous orientation. Quaternion signs
are immaterial; no orientation fitting from deposited spheres is performed.

The point is the transformed **original STL origin**, not its COM. With original
vertex v, uniform scale s, original-frame scaled target COM c and accepted world
COM C, the writer emits P=C−R c. The viewer therefore reconstructs

```text
v_world = P + R (s v) = C + R (s v − c).
```

Scale acts first, then rotation, then translation. Source frames also preserve
geometry-level translation, rotation and scaling before deposition. The COM and
physical properties remain target-derived; source sphere overlaps are irrelevant
to this mass distribution.

`run.conf` records the unchanged original STL in a `# source_stl` comment, as a
quoted path relative to the run directory. Source checksum, declared units, scale
and physical properties are retained in the same record; no sidecar is written. Open
that source STL alongside the pose dataset and click **Apply** on both readers.
Select the pose dataset's **3D Glyphs** representation. Clear the Properties search
box and enable advanced properties with the gear icon. Choose **Pipeline Connection**
for Glyph Type, then select the original STL in the **Input** dropdown immediately
below Glyph Type (initially `(none)`). No dialog opens when selecting Pipeline
Connection. Enable orientation, select `orientation`, set
Orientation Mode to **Quaternion**, and enable scaling using `scale`, Magnitude,
Scale Factor **1**. Hide the standalone STL actor; color instances by `particle_id`.
Do not center or independently rescale the source mesh. The mapper instances one
source mesh; neither file contains a per-particle triangulated copy.

The four components are a point-associated FIELD array, not a `VECTORS` record
(which permits only three). ParaView's quaternion mode uses WXYZ; its OpenGL
mapper scales rotation-matrix columns before adding point translation. See the
[representation definition](https://github.com/Kitware/ParaView/blob/v6.0.1/Remoting/Views/Resources/views_and_representations.xml)
and [VTK mapper implementation](https://github.com/Kitware/VTK/blob/v9.5.2/Rendering/OpenGL2/vtkOpenGLGlyph3DMapper.cxx).

## Source provenance and limitations

New STL-generated particle snapshots use `SDBD_MULTISPHERE 3`, extending version 2
with `source_stl`, `source_origin` and `source_orientation_wxyz`. File I/O preserves
the source frame and writes relocatable references where paths share a root.
Version-1/2 input remains readable for deposition. Older snapshots lack enough
information to establish the original mesh frame, so pose export explicitly
requires regeneration from the original STL. A missing or checksum-mismatched
source is rejected; no mesh location or orientation is guessed. Collision-only
replay does not require the original STL to remain present.

The current engine repeats one geometry per run; each pose dataset consequently
references one original STL. There is no mixed-geometry simulation mode or mesh
registry. Separate geometries/runs have separate output directories/datasets.

STL instances remain continuous rigid bodies across a periodic boundary, whereas
constituent coordinates in `deposition.vtk` are individually wrapped. Equivalent
periodic images can therefore look different near a cell edge. Original STL
surfaces can intersect despite valid multisphere separation when the collision
approximation omits parts of the target. Export does not repair or change physics.
The original contact-position correction also permits small nonrigid numerical
drift of individual spheres; a rigid pose preserves the accepted COM/rotation,
not those per-sphere corrections. Existing tolerances and trajectories are retained.

The WXYZ glyph convention was verified with ParaView 6.0.1. Coordinates and
properties retain solver precision in text; rendered transforms are lower precision. Off-origin meshes can require a wider camera
and clipping range than point bounds suggest. Adjust the view rather than
recentering the mesh or modifying simulation coordinates.
