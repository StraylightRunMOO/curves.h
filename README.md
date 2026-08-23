# curves.h

A single-file, header-only C99/C11 library of **space-filling curves, Gray
codes, and spatial orderings** — locality-preserving maps between 1D indices
and 2D/3D/n-D integer grids.

Drop `include/curves.h` into your project and `#include "curves.h"`. No
allocation, no I/O, no floating-point in the discrete maps, no dependencies
beyond `<stdint.h>` and `<stddef.h>`. Every function is `static inline`
(overridable via `CURVE_INLINE`) and C++-safe.

```c
uint32_t d = curve_hilbert_xy2d(3, 5, 8);   /* 52 on a 256×256 grid */
uint32_t x, y;
curve_hilbert_d2xy(d, 8, &x, &y);           /* x=3, y=5 */
```

## Features

| Family | Functions | Notes |
|---|---|---|
| Gray code | `curve_gray_encode`, `curve_gray_decode` | binary-reflected, 64-bit |
| Morton / Z-order 2D | `curve_morton2`, `curve_morton2_decode` | 16 bits/axis → `uint32_t` |
| Morton / Z-order 2D, wide | `curve_morton2_64`, `curve_morton2_64_decode` | 32 bits/axis → `uint64_t` |
| Morton / Z-order 3D | `curve_morton3`, `curve_morton3_decode` | 21 bits/axis → `uint64_t` |
| Hilbert 2D | `curve_hilbert_xy2d`, `curve_hilbert_d2xy` | constant-time, `bits` 1..16 |
| Hilbert 2D, wide | `curve_hilbert2_64_xy2d`, `curve_hilbert2_64_d2xy` | `uint64_t` index, `bits` 1..32 |
| Hilbert 3D | `curve_hilbert3_xyz2d`, `curve_hilbert3_d2xyz` | A4 lookup tables, `bits` 1..10 |
| Hilbert n-D | `curve_hilbertnd_xy2d`, `curve_hilbertnd_d2xy` | `coords[]` arrays, `n_dims*bits ≤ 64` |
| Moore | `curve_moore_xy2d`, `curve_moore_d2xy` | cyclic Hilbert variant, `bits` 1..16 |
| Sierpiński | `curve_sierpinski_xy2d` | ordering key (non-dense, see below), `bits` 1..31 |
| Windowed | `curve_hilbert_window`, `curve_z_window` | `double` coords through `{xmin,ymin,xmax,ymax}` onto a 16-bit grid |
| Gilbert 2D/3D | `curve_gilbert2`, `curve_gilbert3` | visitor-callback traversal of **arbitrary** rectangles/boxes |
| Batch | `curve_*_batch`, `curve_simd_isa` | SoA arrays; NEON / AVX2 / AVX-512 when compiled in |

Conventions shared by the whole API:

- `bits` selects the grid resolution: side length `2^bits`, index range
  `2^(dims*bits)`. Out-of-range `bits` values are **clamped, never asserted**.
- Coordinates on Hilbert / Moore / Sierpiński / n-D maps are **masked** onto
  that grid (`x & (2^bits - 1)`), so extra high bits wrap rather than
  pollute the index.
- `xy2d` maps coordinates → index; `d2xy` maps index → coordinates.
- Decode output pointers are `restrict` (they must not alias).
- Bit interleaving uses BMI2 `PDEP`/`PEXT` on the scalar path when compiled
  with `-mbmi2`. Batch kernels use SIMD SWAR (NEON on aarch64, AVX2/AVX-512
  on x86). Define `CURVE_NO_SIMD` to force scalar batch loops.

## Quick start

```c
#include "curves.h"
#include <stdio.h>

int main(void)
{
    uint32_t d = curve_hilbert_xy2d(3, 5, 8);
    uint32_t x, y;
    curve_hilbert_d2xy(d, 8, &x, &y);

    /* Many points at once (SoA, no aliasing): */
    uint32_t xs[4] = {0, 1, 2, 3}, ys[4] = {0, 0, 1, 1}, zs[4];
    curve_morton2_batch(xs, ys, zs, 4);

    printf("%u %s\n", d, curve_simd_isa());
    return 0;
}
```

```bash
gcc -std=c11 -O2 -I include main.c -o main
# x86: add -mbmi2 -mavx2  (or -mavx512f) for the fast paths
# ARM: NEON batch kernels are used automatically on aarch64
```

## Notes on specific curves

- **Moore** is the closed-loop variant of Hilbert: the last cell is adjacent
  to the first. Order 1 equals order-1 Hilbert; for `bits ≥ 2` each quadrant
  is a rotated Hilbert sub-curve.
- **Sierpiński** is an *ordering key*, not a dense index: the curve is closed
  and its two triangular halves share the main diagonal, so keys live in
  `[0, 2^(2*bits+1))` with half the space unused. Keys are injective over
  the grid; sorting points by key yields Sierpiński order (useful e.g. for
  TSP-style heuristics, Bartholdi & Platzman 1989). No `d2xy` inverse.
- **Hilbert n-D** (Skilling's formulation of the Butz algorithm) takes a
  coordinate array: `curve_hilbertnd_xy2d(coords, n_dims, bits)`. `n_dims < 1`
  is a no-op (returns 0 / writes nothing). For `n_dims = 2` it matches
  `curve_hilbert_xy2d`; for `n_dims = 3` the orientation convention differs
  from `curve_hilbert3_xyz2d` (both are valid Hilbert curves).
- **Gilbert** fills arbitrary rectangles/boxes via a visitor callback:
  `curve_gilbert2(w, h, cb, user)` calls `cb(x, y, user)` once per cell in
  curve order. 2D paths are always edge-adjacent; 3D paths are face-adjacent
  when all box dimensions are even (with odd dimensions, joins between
  sub-blocks may step diagonally — as in the reference implementation).
- **Windowed** helpers map floating-point coordinates inside
  `window = {xmin, ymin, xmax, ymax}` onto the 16-bit integer grid. Zero-size
  axes map to 0. NaN and −Inf clamp to 0; +Inf clamps to `0xFFFF`.

## Batch API

Every batch function takes structure-of-arrays buffers that **must not
alias**, plus a count. Unaligned pointers are fine.

```c
void curve_morton2_batch(const uint32_t *x, const uint32_t *y,
                         uint32_t *out, size_t n);
void curve_morton2_decode_batch(const uint32_t *d, uint32_t *x, uint32_t *y,
                                size_t n);
void curve_morton2_64_batch(...);     /* uint32 x,y → uint64 */
void curve_morton3_batch(...);        /* uint32 x,y,z → uint64 */
void curve_gray_encode_batch(...);
void curve_hilbert_xy2d_batch(const uint32_t *x, const uint32_t *y,
                              int bits, uint32_t *out, size_t n);
```

`curve_simd_isa()` returns `"neon"`, `"avx2"`, `"avx2+bmi2"`, `"avx512"`,
`"avx512+bmi2"`, `"bmi2"`, or `"scalar"`.

32-bit Morton 2D batch is vectorized on NEON, AVX2, and AVX-512. 64-bit
Morton 2D/3D batch uses AVX2/AVX-512; on NEON those fall back to a scalar
loop (2-wide SWAR lost to the compiler).

## Building the tests, benchmarks, and examples

Direct compile (no CMake needed):

```bash
gcc -std=c11 -pedantic -Wall -Wextra -O2 -I include tests/test_curves.c -lm -o test_curves
./test_curves
```

With CMake:

```bash
cmake -S . -B build -DCURVES_BUILD_BENCHMARKS=ON
cmake --build build
ctest --test-dir build --output-on-failure
./build/benchmarks/bench_main          # throughput numbers
./build/examples/render_pgm /tmp/out   # morton/hilbert/moore/sierpinski/gilbert .pgm
./build/examples/locality_demo         # row-major vs Morton vs Hilbert vs Moore
```

| Option | Default | |
|---|---|---|
| `CURVES_BUILD_TESTS` | ON | C suite + C++11 smoke test |
| `CURVES_BUILD_EXAMPLES` | ON | `render_pgm`, `locality_demo` |
| `CURVES_BUILD_BENCHMARKS` | OFF | `bench_main` (`-O3 -march=native`) |
| `CURVES_ENABLE_SANITIZERS` | OFF | ASan + UBSan, no recover |

## CMake integration

```bash
cmake --install build --prefix /your/prefix
```

```cmake
find_package(curves REQUIRED)
target_link_libraries(your_target PRIVATE curves::curves)
```

or `add_subdirectory(curves.h)` and link `curves::curves` the same way.

## Provenance

- 2D Hilbert (constant-time) and windowed maps:
  [rawrunprotected/hilbert_curves](https://github.com/rawrunprotected/hilbert_curves)
  and [tidwall/curve](https://github.com/tidwall/curve) (both public domain).
- 3D Hilbert A4 lookup tables: same source.
- Gilbert curves: [jakubcerveny/gilbert](https://github.com/jakubcerveny/gilbert)
  (BSD-2-Clause) via skookum/spatial.
- n-D Hilbert: A. R. Butz, *Alternative Algorithm for Hilbert's Space-Filling
  Curve* (IEEE Trans. Comp., 1971); Doug Moore's Rice University
  implementation; formulation after J. Skilling, *Programming the Hilbert
  curve* (2004).
- Moore curve: E. H. Moore (1900); quadrant structure derived from the
  curve's Lindenmayer system.
- Sierpiński ordering: J. J. Bartholdi III & L. K. Platzman, *Spacefilling
  curves and the planar traveling salesman problem* (JACM, 1989).

## License

MIT — see [LICENSE](LICENSE).
