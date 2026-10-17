# Cube: original surface-lattice representation

From the repository root:

```sh
build/deposit run examples/cube/config.conf
```

This small example deposits **12 bodies**, seed **84328**, GSL **mt19937**, in a
**1440 × 1440 × 1440 mm** observation cell with periodic sides and an open top.
Each nominal **180 mm** cube uses the original perturbed surface lattice at
resolution **3**: **26 overlapping spheres**, rather than a solid exact cube.
The original COM approximation and generator random draws are retained.

Outputs: `results/cube/run.conf` and `results/cube/deposition.vtk`.
Use the main README's [ParaView sphere settings](../../README.md#build-run-and-view).
The example favors a quick run rather than a tall pile. Change `count` and `box`
in a copy of the config for a larger study, choosing a fresh output directory.
