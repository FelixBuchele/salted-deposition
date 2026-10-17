# Cube deposition

This example deposits 12 cube-shaped particles using the sphere representation
from Nikola Topic's original implementation. It provides a simple starting point
and a geometry for comparisons with the original code.

## Run

After [building SALTED](../../README.md#build-salted), run from the project folder:

```bash
build/deposit run examples/cube/config.conf
```

The settings in [config.conf](config.conf) specify a nominal cube edge of 180 mm,
resolution 3 and a 1440 × 1440 × 1440 mm domain. Each particle has 26 spheres on
its surface. The random seed is 84328 and the generator is mt19937. Horizontal
boundaries are periodic; the maximum packing height is 1440 mm.

## Results and customization

Results are written to `results/cube/`:

- `deposition.vtk`: the deposited spheres, grouped by particle identifier.
- `run.conf`: the simulation settings and information about the completed run.

Use the VTK file to inspect particle orientations and contacts. See the
[export guide](../../docs/exports.md) for visualization settings.

Change `count` and `box` for a larger deposit, or `size` and `resolution` to change
the particle dimensions and sphere discretization. Choose a new output folder
with `--output` when repeating a run to preserve previous results.

The cube retains the original small geometric perturbations and center-of-mass
approximation. Its sphere surface approximates a cube; it does not have exact
flat faces. See the [particle reference](../../docs/particles.md#generators) for details.
