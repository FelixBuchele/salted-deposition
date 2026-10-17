<div align="center">
<h1><img src="docs/salted_logo.png" alt="SALTED — Simulation of Granular Packings via Sequential Particle Deposition" width="800"></h1>
<img alt="C++20" src="https://img.shields.io/badge/C%2B%2B-20-blue">
<a href="LICENSE"><img alt="SALTED source: MIT" src="https://img.shields.io/badge/source-MIT-blue"></a>
</div>

---

SALTED constructs static deposits of rigid multisphere particles by sequential
steepest-descent ballistic deposition [1]. Each body settles through gravity-driven
contact and rotation events before becoming fixed. The substrate is z=0 and
horizontal boundaries are periodic; this is geometrical deposition rather than
time-resolved force-based DEM.

<p align="center">
<img src="docs/meshes_03.png" alt="Steepest-descent deposition of 150 cow2 particles, rendered using their original STL surfaces" width="720">
</p>

*Figure 1. Steepest-descent deposit of 150 cow2 particles, rendered using original
STL instances. Each particle uses 2,189 spheres for collision calculations.
Cow model: Viewpoint Animation Engineering, Sun Microsystems;
[Princeton Suggestive Contour Gallery](https://gfx.cs.princeton.edu/proj/sugcon/models/).
The source mesh is not distributed with SALTED and retains its creators' rights.*

## Features

- Built-in sphere, spiky, original surface-lattice cube and overlapping-sphere stick.
- Manual particle files with explicit physical centers of mass (COM).
- Bundled, optional GEMSS generation from binary STL [2,3]. Physical COM, volume
  and inertia come from the target voxel geometry (`compute_physics=2`),
  independently of overlapping collision spheres.
- Runtime configuration, reproducible seeds, replay metadata and optional Python.
- Compact VTK collision geometry and original-STL instances for ParaView.

## Build, run and view

Requirements: Linux/glibc, C++20 compiler, CMake ≥3.20, GSL ≥2.0 development
files and a build tool. From the repository root:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 2
build/deposit run examples/cube/config.conf
```

The run writes **`results/cube/deposition.vtk`** for visualization and **`run.conf`**
with effective settings and provenance. See [output and visualization](#output-and-visualization)
for the exported geometry and its use in tools such as ParaView.

The default build needs neither Python nor Eigen and runs the prepared rock.
Existing output destinations are refused. Repeat with a fresh destination:

```sh
build/deposit run examples/cube/config.conf --output results/cube-repeat
```

## Configuration

Examples use `key = value` files. Edit `count`, `seed`, `particle` and
`box = Lx Ly height`; lengths use the declared `units`. The top is open: height
limits the accepted surface height rather than adding a lid.

```sh
build/deposit validate --config examples/cube/config.conf
build/deposit run --help
```

Validation writes nothing and does not simulate or reconstruct. Precedence is
**defaults < config < --overrides FILE < --output**. Config paths are file-relative;
`--output` is cwd-relative. See the [configuration reference](docs/configuration.md)
and `--help`, `--version`, `--print-example` for command discovery.

## Python interface (optional)

Python ≥3.10; no runtime packages. Install with `python3 -m pip install ./python`
(editable: `python3 -m pip install -e ./python`). Build C++ separately:

```python
from salted_deposition import Simulator

sim = Simulator("build/deposit")
result = sim.run("examples/cube/config.conf", output="results/python-cube",
                 count=4, seed=12345)
print(result.vtk_path)
```

The API exposes `Simulator`, `run()`, `generate_particle()` and immutable
`SimulationResult`. Common overrides are `count`, `seed`, `shape` or `particle`,
and `export_stl_poses`; None retains config values. VTK is automatic. Advanced
settings remain in configs. See [API details](docs/configuration.md#python-api).

## STL generation and reuse

For reconstruction, provide Eigen ≥3.4 and enable the bundled GEMSS backend:

```sh
cmake -S . -B build-gemss -DCMAKE_BUILD_TYPE=Release -DSDBD_ENABLE_GEMSS=ON
cmake --build build-gemss --parallel 2
build-gemss/deposit generate examples/rock/generation.conf
printf 'particle = rock-regenerated.particle\n' > results/rock-input.conf
build/deposit run examples/rock/config.conf --overrides results/rock-input.conf --output results/rock-regenerated
```

The overlay selects the generated file relative to `results/rock-input.conf`.
Normal rock runs reuse the included `rock.particle`. CMake finds installed Eigen;
`SDBD_EIGEN_INCLUDE_DIRECTORY` may instead name a local header root containing
`Eigen/Core`. Builds never download dependencies. Python generation and reuse:

```python
from salted_deposition import Simulator

sim = Simulator("build-gemss/deposit")
particle = sim.generate_particle("examples/rock/example_mesh.stl",
                                 output="results/python-rock.particle",
                                 config="examples/rock/generation.conf")
for seed in (84328, 84329):
    result = sim.run("examples/rock/config.conf", particle=particle,
                     output=f"results/python-rock-{seed}", count=6, seed=seed,
                     export_stl_poses=True)
    print(result.vtk_path, result.stl_poses_path)
```

STL units must be declared. [Particle formats](docs/particles.md) cover manual
input; [rock notes](examples/rock/README.md) give the Dice **0.9** generation settings.

## Examples

| Example | Purpose | Command after `build/deposit` |
|---|---|---|
| [Cube](examples/cube/README.md) | 12 bodies; original 26-sphere surface lattice | `run examples/cube/config.conf` |
| [Stick](examples/stick/README.md) | 12 bodies; seven overlapping spheres each | `run examples/stick/config.conf` |
| [Rock](examples/rock/README.md) | 24 bodies; reusable STL geometry and surface poses | `run examples/rock/config.conf` |

These are small demonstrations, not large pile showcases. They pin GSL mt19937
and seeds. Notes give dimensions and expected files. The rock mesh belongs to
Felix Buchele and is included under MIT; no large outputs are distributed.

## Output and visualization

| File | Contents |
|---|---|
| `run.conf` | Replay settings, source geometry and units, status, final physical COMs and redrop counts |
| `deposition.vtk` | Final world-space sphere centers, `radius`, `particle_id` and local `sphere_id` |
| `stl_poses.vtk` | STL-derived bodies: transformed original mesh origins, `particle_id`, WXYZ `orientation` quaternions and uniform `scale` |

Both VTK files are ASCII PolyData. `deposition.vtk` stores one point per constituent
sphere, rather than triangulated sphere surfaces. Visualization tools can render
spheres at these points using their radii to reconstruct the collision geometry.
Constituents sharing a `particle_id` belong to the same rigid body; `sphere_id`
identifies each constituent within it. This view supports inspection of packing
structure and contacts, with particles colored by identity.

`stl_poses.vtk` stores one point per deposited body and requires the **unchanged
original STL** referenced in `run.conf`. Each point is the transformed mesh origin,
**not the physical COM**. Its normalized, right-handed WXYZ quaternion and scale
allow a tool such as ParaView to scale, rotate and place instances of that mesh,
reconstructing the final surface assembly without duplicating mesh data in the
export. Coordinates use the simulation's declared length unit; `scale` converts
original STL coordinates to that unit, including any specified enlargement.
The surface view shows particle shapes and orientations; the sphere view shows
the geometry actually used for collision calculations. Approximation differences
can therefore produce visible STL intersections despite separated collision spheres.

`stl_poses = false` disables poses; `legacy_output = true` adds original numerical
text files. Replay with
`build/deposit run results/cube/run.conf --output results/cube-replay`.
[Export conventions](docs/exports.md) give coordinate transforms, precision and
optional ParaView setup instructions.

## Limitations

One geometry is repeated per run. Horizontal box lengths must be positive integers
in the declared unit. Absolute tolerances and original stopping rules remain:
scale changes can change trajectories, some shapes/seeds stall, and redrop does
not guarantee termination. Exact comparisons depend on long-double precision and
compiler/libm. Windows/macOS are not currently supported.

Voxelization can omit features; mesh checks do not certify absence of
self-intersections. STL surfaces can intersect where collision spheres omit
geometry. Periodic wrapping and legacy contact corrections can make rigid poses
differ slightly from individual sphere positions. ParaView arithmetic has lower
precision than the solver. See [export limitations](docs/exports.md#source-provenance-and-limitations).

## License, authors and references

SALTED source retains the [MIT license](LICENSE), under maintainer-confirmed
original-source permission. Original code: **Nikola Topic** (2015); methodology:
**Nikola Topic and Thorsten Pöschel** [1]; modernization: **Felix Buchele**.
Original vector code credits **Thomas Schwager**. GEMSS author **Arash Moradian**
and contributors **Felix Buchele, Patric Müller and Thorsten Pöschel** retain
multisphere-generation attribution.

GSL is GPL-3.0-or-later. Bundled GEMSS is pinned to
`dc0f6de101b04545aa9679a13e692f60dfae2689`; MIT does not replace its component
MPL/LGPL/MIT/zlib notices. Eigen remains external. Read
[bundled notices](vendor/GEMSS/NOTICE.md) and [distribution terms](docs/licenses.md).
Linked binaries are not MIT-only and require applicable GPL/LGPL/MPL compliance.
This source distribution includes example geometry, not compiled executables.

1. N. Topic and T. Pöschel. *Steepest descent ballistic deposition of complex
   shaped particles*. Journal of Computational Physics **308**, 421–437 (2016).
   [DOI](https://doi.org/10.1016/j.jcp.2015.12.052).
2. A. Moradian, F. Buchele and T. Pöschel. *GEMSS: A C++ Library for Multi-Sphere
   Modeling in DEM Simulations*. arXiv:2608.15463 (2026).
   [Verified preprint](https://arxiv.org/abs/2608.15463). See also the
   [GEMSS repository citation](https://github.com/aqa-arash/GEMSS#citation).
3. F. Buchele, T. Pöschel and P. Müller. *Multi-sphere shape generator for DEM
   simulations of complex-shaped particles*. Powder Technology **480**, 122635
   (2026). [DOI](https://doi.org/10.1016/j.powtec.2026.122635).
