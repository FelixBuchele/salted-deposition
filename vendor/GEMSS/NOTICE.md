# Bundled GEMSS source

Source: https://github.com/aqa-arash/GEMSS

Revision: `dc0f6de101b04545aa9679a13e692f60dfae2689`.

The 149 source files are the complete local include closure of SALTED's GEMSS
adapter with Eigen 3.4.0 on GNU/Linux, including locally referenced conditional
includes. Paths below this directory preserve upstream `src/GEMSS` layout. The
upstream application, examples, unrelated libigl modules and build system are
omitted. No copied upstream source file has been modified. `manifest.sha256`
records the upstream files and license texts byte for byte. Eigen is external.

| Component | Authorship and terms |
|---|---|
| GEMSS headers | Arash Moradian, Felix Buchele, Institute for Multiscale Simulation (FAU), 2025–2026; MIT, `LICENSE` |
| libigl subset | Individual copyright notices retained in each file; MPL-2.0, `thirdparty/igl/LICENSE.MPL2` |
| `FastWindingNumberForSoups.h` | Side Effects Software Inc. (2018); MIT notices embedded in the file |
| `raytri.c` | Tomas Moller (2000), with libigl alterations; retained upstream public-domain statement |
| `edt.hpp` | William Silversmith, Seung Lab; LGPLv3, with governing `thirdparty/COPYING.LESSER` and `thirdparty/COPYING` |
| `threadpool.h` | Jakob Progsch and Václav Zeman (2012), with William Silversmith alterations; zlib-style notice embedded in the file |

The bundled tree is not uniformly MIT. Source availability under the MPL and LGPL
is provided by these unmodified source files; their notices must accompany
redistribution. Binary distribution requires compliance with applicable component
terms in addition to GSL's GPL obligations. No linked binaries are included.

The full upstream libigl tree contains additional copyleft/optional modules;
these are not included. This subset is sufficient for SALTED's configured
header-only use, not a general libigl installation. Builds use
`EIGEN_MPL2_ONLY`, no OpenMP linkage and no native/fast-math flags.
