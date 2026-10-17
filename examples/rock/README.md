# Rock: reusable STL-derived particle

Felix Buchele owns the unchanged `example_mesh.stl` and authorized MIT distribution
of this mesh and `rock.particle`. Coordinates are **mm**, uniform scale **1**.
The binary STL contains **96,250 triangles**, dimensions approximately
**9.252 × 9.998 × 7.275 mm**; its hash is recorded in the license notes.

From the repository root:

```sh
build/deposit run examples/rock/config.conf
```

The prepared file needs no GEMSS build. It deposits **24 bodies**, seed **84328**,
GSL **mt19937**, in a **60 × 60 × 60 mm** observation cell, retaining the original
phi/redrop defaults. Outputs are `results/rock/run.conf`, `deposition.vtk` and
`stl_poses.vtk`, all in the same directory.

## Regenerate and reuse

Build with GEMSS as described in the main README, then:

```sh
build-gemss/deposit validate --config examples/rock/generation.conf
build-gemss/deposit generate examples/rock/generation.conf
printf 'particle = rock-regenerated.particle\n' > results/rock-input.conf
build/deposit run examples/rock/config.conf --overrides results/rock-input.conf --output results/rock-regenerated
```

The generated file is `results/rock-regenerated.particle`; the overlay selects
it relative to `results/rock-input.conf`. Choose new output paths for repeated runs.
These commands leave the included prepared file unchanged.

Settings: voxel resolution **24**, cap **50 spheres**, requested Dice **0.9**,
minimum radius **1 voxel**, relative center distance **0.6**. Other backend
controls retain their defaults; confinement and pruning are disabled. The
prepared reconstruction contains **19 spheres**, achieved voxel Dice **0.900792**,
and spacing **0.303134 mm**. It is an economical approximation with omitted fine
features, not an exact surface reconstruction.

Physical properties use the original target voxel grid (`compute_physics=2`),
independently of collision-sphere overlaps. Stored COM is approximately
**(-0.044523, -0.150169, -0.383261) mm**, and target volume **247.242584 mm³**.
The source origin and frame are preserved.

## View the deposition

Use [ParaView 3D Glyphs](../../docs/exports.md) for collision spheres or original
STL instances. For the latter, open `example_mesh.stl` and `stl_poses.vtk`, Apply
both, select the poses, clear the Properties search and enable advanced properties
(gear). Choose **Glyph Type → Pipeline Connection**, then select the STL in the
**Input** dropdown below it. Selecting Pipeline Connection opens no dialog.
Use WXYZ Quaternion orientation and the `scale` array with factor **1**.
Do not recenter the STL. Contacts use the sphere union; visible STL intersections
need not mean sphere penetration.
