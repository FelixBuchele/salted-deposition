# Rock deposition from STL

This example uses an irregular rock defined by [example_mesh.stl](example_mesh.stl).
The included [rock.particle](rock.particle) already contains its 19-sphere
representation and physical properties, so the simulation can run without
repeating GEMSS generation.

## Run

After [building SALTED](../../README.md#build-salted), run from the project folder:

```bash
build/deposit run examples/rock/config.conf
```

The settings in [config.conf](config.conf) deposit 24 rocks in a 60 × 60 × 60 mm
domain with periodic horizontal boundaries. The seed is 84328 and the generator
is mt19937. The mesh coordinates are millimetres and the scale factor is 1.

## Results

Results are written to `results/rock/`:

| File | Use |
|---|---|
| `deposition.vtk` | Display the collision spheres and inspect contacts. |
| `stl_poses.vtk` | Combine with the unchanged `example_mesh.stl` to display the deposited rock surfaces. |
| `run.conf` | Inspect or reuse the simulation settings and recorded physical centers of mass. |

The sphere view shows the geometry used by the solver. The STL view shows the
original rock surface, including details omitted by the approximation. See the
[export guide](../../docs/exports.md) for formats and ParaView setup.

## Regenerate or use another mesh

Enable GEMSS using the [build instructions](../../README.md#build-salted), then run:

```bash
build/deposit validate --config examples/rock/generation.conf
build/deposit generate examples/rock/generation.conf
```

The [generation configuration](generation.conf) uses voxel resolution 24, at most
50 spheres and requested Dice similarity 0.9. The included particle file records
19 spheres and achieved Dice 0.900792. Its remaining generation parameters are
stored in the same file; the generation configuration sets them for regeneration.

Generation writes `results/rock-regenerated.particle`, leaving the included
`rock.particle` unchanged. Use the complete `rock-study.conf` example in the main
README's [custom-geometry workflow](../../README.md#use-your-own-particle-geometry)
to deposit the regenerated particle. For another mesh, copy `generation.conf`
and change its `stl`, `units` and `output` settings. Choose new output names for
repeated generation or simulation runs.

Physical volume, center of mass and inertia come from the voxelized original
mesh, independently of collision-sphere overlaps. These properties and the source
mesh reference are preserved in the particle file.

The mesh and generated particle definition are owned by Felix Buchele and
included under MIT. See the [license notes](../../docs/licenses.md) for attribution
and the original mesh checksum.
