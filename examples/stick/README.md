# Stick: elongated overlapping-sphere particle

From the repository root:

```sh
build/deposit run examples/stick/config.conf
```

This small example deposits **12 bodies**, seed **4**, GSL **mt19937**, in a
**960 × 960 × 960 mm** observation cell with periodic sides and an open top.
Each stick has total length **180 mm**, including end caps, and **seven spheres
of radius 30 mm** along its body x-axis. Centers are **20 mm** apart, so adjacent
constituents overlap. Uniform union symmetry places the physical COM at the origin.
The normal solver controls each body's placement, rotation and settlement.

Outputs: `results/stick/run.conf` and `results/stick/deposition.vtk`.
Use the main README's [ParaView sphere settings](../../README.md#build-run-and-view).
This is a compact overlap/orientation demonstration, not a dense pile showcase.
Adjust `length`, `radius` and `spheres` together to preserve positive internal overlap.
