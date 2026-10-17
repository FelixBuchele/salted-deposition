# Stick deposition

This example deposits 12 elongated particles, each constructed from seven
overlapping spheres. It demonstrates orientation-dependent settling with a
customizable particle shape.

## Run

After [building SALTED](../../README.md#build-salted), run from the project folder:

```bash
build/deposit run examples/stick/config.conf
```

The settings in [config.conf](config.conf) define sticks 180 mm long, including
the end caps, with sphere radii of 30 mm and center spacing of 20 mm. The domain
is 960 × 960 × 960 mm with periodic horizontal boundaries. The random seed is 4
and the generator is mt19937.

At uniform density, the stick's symmetry places its physical center of mass at
its geometric center. Internal sphere overlaps are part of the same rigid object.

## Results and customization

Results are written to `results/stick/`:

- `deposition.vtk`: the final sphere geometry, including internal overlaps.
- `run.conf`: the simulation settings and information about the completed run.

Use the VTK file to inspect stick orientations and contacts. See the
[export guide](../../docs/exports.md) for visualization settings.

Change `length`, `radius` and `spheres` to define another stick. The center spacing
`(length - 2*radius)/(spheres - 1)` must be positive and smaller than `2*radius`,
so adjacent spheres overlap. Change `count` and `box` to control the deposit.
Use `--output` to select a new folder for each run.
