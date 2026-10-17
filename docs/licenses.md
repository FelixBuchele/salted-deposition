# Source and third-party terms

The included MIT notice preserves Nikola Topic's original implementation
copyright and Felix Buchele's modernization credit. Public source redistribution
under these terms was authorized by the maintainer. The supplied rock
`examples/rock/example_mesh.stl` is owned by Felix Buchele; he explicitly authorized
its inclusion and the generated `rock.particle` under the same MIT terms on
17 October 2026. Mesh coordinates are millimetres. The original mesh bytes are
unchanged: SHA-256 `928773dbf822563f71e5681f68e6bfd7051a853636536cb7ca37f5dbe4aa983e`.
No other user-supplied mesh, PDF, archive or generated derivative is distributed.

Required GSL is [GPL-3.0-or-later](https://www.gnu.org/software/gsl/).
The optional bundled [GEMSS subset](../vendor/GEMSS/NOTICE.md) is version-pinned,
unmodified and accompanied by its original notices and source hashes. It includes:

| Component | Governing notices |
|---|---|
| GEMSS | MIT, `vendor/GEMSS/LICENSE` |
| Used libigl headers | MPL-2.0, individual notices and `thirdparty/igl/LICENSE.MPL2` |
| Embedded SideFX winding code | MIT notice in `FastWindingNumberForSoups.h` |
| EDT by William Silversmith | LGPLv3; `thirdparty/COPYING.LESSER` and `COPYING` |
| Thread pool | Embedded Progsch/Zeman zlib-style notice and Silversmith alteration notice |
| `raytri.c` | Retained Moller/libigl public-domain statement |

Eigen is external (MPL/permissive); builds enforce `EIGEN_MPL2_ONLY`. Unrelated
libigl copyleft modules and upstream applications are not bundled. Neither GEMSS's
MIT notice nor SALTED's MIT notice replaces component terms. Distribution of
GSL-linked executables requires GPL compliance, and header-instantiated GEMSS
components bring their corresponding LGPL/MPL obligations. Preserve applicable
notices and provide corresponding source as required; do not label linked binaries
MIT-only. No compiled executables are included. Python packaging adds no runtime dependency.
