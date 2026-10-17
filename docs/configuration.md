# Configuration reference

The public executable is `deposit`. Use `run CONFIG`, `generate CONFIG` or
`validate --config CONFIG`; positional `validate CONFIG` remains accepted.
`COMMAND --help` shows operation-specific help. `--print-example` writes a
starter configuration to stdout. Successful commands return 0, handled errors 1.

## Parsing and precedence

One `key = value` per line; `#` begins a comment outside double quotes.
Unknown/duplicate keys, trailing numeric text and misplaced settings fail.
Booleans are exactly `true` or `false`. `format = SALTED_2` is optional for
handwritten files and is written in saved run records. Particle files use a
different format with fields in a fixed order; see [particle formats](particles.md).

Precedence: **defaults < CONFIG < --overrides FILE < --output PATH**.
Each path resolves beside the file that declares it; `--output` resolves from
the current working directory. Paths passed to Python methods also use the current
directory and become absolute before the command runs. Original configurations
are not rewritten; C++ reads and validates all settings.

All operations accept one `--overrides FILE`; `run` and `generate` accept
`--output PATH`. An overlay replacing `particle` drops a recorded old geometry
checksum; an explicitly supplied new checksum is still verified. Other keys
remain untouched, including incompatible geometry-specific settings. Do not use
a cube `size` override with a file particle or a stick. `--config FILE` remains
accepted wherever a positional configuration is accepted; prefer the forms above.

Output destinations must be new. `run` and `generate` require `output` unless
`--output` supplies it. Validation does not need output, creates no files,
consumes no random draws and does not simulate or reconstruct. Mesh validation
can load/check topology and resource estimates but cannot guarantee successful
volumetric reconstruction.

## Simulation

Only the destination is required; other keys retain these existing defaults:

| Key | Meaning / default |
|---|---|
| `particle` | `sphere`, `spiky`, `cube`, `stick` or prepared file; default `sphere` |
| `count` | Rigid-particle count, 24; not constituent-sphere count |
| `seed` | GSL seed, 84328 |
| `rng` | Explicit GSL name; otherwise `GSL_RNG_TYPE`, otherwise mt19937 |
| `box` | x/y periods and maximum accepted surface height, 720 720 6680 |
| `units` | `mm` (default), `cm`, `m`, `legacy`; labels all run lengths |
| `size` | Built-in nominal outer size; default sphere diameter 60, cube edge 60·resolution, spiky extent 60·(2·resolution−1) |
| `resolution` | Original built-in discretization: cube 3, otherwise 1; cube requires ≥2 |
| `length`, `radius`, `spheres` | Stick total length, constituent radius and count: 180, 30, 7 |
| `phi` | Small rotation increment used in contact-release tests, 0.01 radians; not a time step |
| `redrop` | Original retry policy, false; no termination guarantee |
| `injection` | `uniform` (default): original full-cell placement; `circle`: centered uniform-area disk |
| `quiet` | Suppress ordinary progress/summary, false |
| `legacy_output` | Add compatibility numerical text files, false |
| `stl_poses` | Export poses when a valid STL source/frame exists, true |
| `particle_checksum` | Saved file-geometry fingerprint; normally leave to replay metadata |
| `output` | Required new result directory |

Original spiky/cube jitter remains, so nominal dimensions are not exact realized
bounds. At resolution 1, spiky degenerates to one sphere; use a higher resolution
to represent its six arms. A resolution-3 cube has 26 surface constituents.
A stick requires 0 < (length−2·radius)/(spheres−1) < 2·radius and spheres ≥2.
Prepared files retain their dimensions and explicit COM. Their declared physical
units must match the run; `legacy` files take the run unit without conversion.
Changing the units label alone does not resize geometry or numerical tolerances.

Numerical collapse means that contact events make too little positional progress
or cannot continue with the current support geometry. With `redrop = true`, the
solver retries with a new horizontal position and orientation; difficult cases
may keep retrying indefinitely.

Sides are periodic, the plane is z=0 and the top is open. Height is an acceptance
limit: exceeding it fails without repositioning. Horizontal lengths must be
positive integers in the declared unit; height may be fractional. The existing
search grid derives spacing `gcd(120, Lx, Ly)` and below-plane margin ≥960.
More than 20 million cells is rejected; commensurate dimensions (e.g. multiples
of 120) avoid unnecessarily fine grids. Fractional horizontal periods and solid
side walls are not supported.

## GEMSS generation

`stl`, declared `units` and an output particle-file path are required.
The backend must be enabled at build time. Preparation is separate from
simulation: reuse its file with `particle = FILE`, even without GEMSS installed.

| Key | Meaning / default |
|---|---|
| `stl` | Binary STL input; source origin and orientation are preserved |
| `units` | Explicit target length unit; also source unit unless overridden |
| `source_units` | Source coordinate unit, default `units` |
| `scale` | Additional uniform enlargement, 1 |
| `voxel_resolution` | Divisions of shortest bounding dimension, 32 |
| `maximum_spheres` | Maximum constituents per rigid body, 64 |
| `accuracy` | Requested Dice reconstruction similarity, 0.95 |
| `minimum_radius_voxels` | Minimum constituent radius in voxels, 2 |
| `minimum_center_distance_relative` | GEMSS relative center-distance parameter, 1.4 |
| `padding`, `search_window`, `persistence` | Backend controls, 2 each |
| `radius_offset_voxels` | Backend radius offset in voxels, 0.5 |
| `confine_mesh`, `prune_isolated_spheres` | Optional backend controls, false |
| `density` | Uniform mass per target length³, 1; choose a consistent mass unit |
| `maximum_voxels` | Resource guard, 5,000,000 |
| `quiet` | Suppress ordinary generation summary, false; accuracy warnings remain |
| `output` | Required new reusable particle file |

Coordinate conversion is `scale × metres(source_units)/metres(units)`; legacy
units cannot be implicitly converted to physical units. Target properties always
use `compute_physics=2`. Voxel-based COM approximates the original STL; no
sphere-volume-weighted fallback exists. More spheres at fixed target resolution
do not redefine the physical COM. Generation reports target COM, spacing and
requested/achieved Dice; missing the target produces a warning, not guaranteed
failure. Inspect geometry before choosing settings for a scientific study.

## Circular injection

`injection = circle` places initial physical COMs in a disk centered at (Lx/2,Ly/2)
with radius R=min(Lx,Ly)/4. For two independent GSL uniforms U1,U2,
r=R sqrt(U1), theta=2 pi U2 and position=(Lx/2+r cos(theta), Ly/2+r sin(theta)).
No radius, bias or center parameters are introduced. Redrop uses the same selected
injection mode. After placement all movement/rotation/contact rules remain the
original ones; particles may settle outside the disk. The mode is written in
`run.conf` with the seed and resolved RNG.

Uniform mode preserves the original expressions, draw order and reference-frame
placement convention, including existing imported-particle COM offsets. Circular
mode intentionally positions the physical COM itself, irrespective of mesh origin.
Both consume two placement uniforms; the changed geometry can cause later
trajectories/random consumption to diverge. No exact compatibility with uniform
trajectories is claimed for circular injection.

## Output and replay

See the main README for output names and replay; [export conventions](exports.md)
define arrays, transforms and precision. Generation accuracy warnings remain
visible even with quiet enabled. Original STL bytes are needed only for pose
export; version-1/2 particle files can be simulated without pose provenance.

## Python API

The optional `salted_deposition` package invokes an already compiled executable.
`Simulator(executable="deposit")` searches PATH; relative executable paths are
anchored at construction. The public methods are:

```text
sim.run(config, *, output, count=None, seed=None, shape=None, particle=None,
        export_stl_poses=None)
sim.generate_particle(mesh, *, output, max_spheres=None, units=None, config=None)
```

`run()` returns immutable `SimulationResult(output_dir, vtk_path, stl_poses_path)`;
the last path is None when unavailable. VTK is automatic. Shape and particle
are mutually exclusive. Replacing geometry retains other configuration fields,
so use a compatible base config. None preserves existing settings.

`generate_particle()` returns the generated file path. Declare units explicitly
or supply a generation config with units. `max_spheres` maps to C++
`maximum_spheres`; advanced controls remain in that config. Python reads no mesh
data and does not calculate physical properties.

Method paths are cwd-relative. Diagnostics stream to the terminal; nonzero exits
raise `subprocess.CalledProcessError`, missing required exports raise
`FileNotFoundError`, and conflicting or invalid API choices raise
`ValueError`/`TypeError`. Existing outputs are refused, matching the CLI.
