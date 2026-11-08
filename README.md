# curves.h

Space-filling curves, Gray codes, and spatial keys. One C99/C11 header.
Drop `include/curves.h` in and `#include "curves.h"`. No heap, no I/O, no
`<math.h>` on the discrete maps, no dependencies besides `<stdint.h>` and
`<stddef.h>`. Every function is `static inline` (overridable via
`CURVE_INLINE`) and C++-safe.

Hilbert will not make two embeddings similar. That is a different library,
and usually a different mistake. This one maps integer grids onto 1-D
indices and keeps nearby cells nearby on the curve. Sockets, databases,
and the reason you still need an exact rerank remain yours.

Version 1.1. The public declarations in the header are the API reference.
We read Skilling so you do not have to. You are welcome.

## Quickstart

The shortest path that does something useful. It does not try to be clever.
It will also not invent a Hilbert curve of its own.

```c
#include "curves.h"
#include <stdio.h>

int main(void)
{
    uint32_t d = curve_hilbert_xy2d(3, 5, 8);   /* 52 on a 256×256 grid */
    uint32_t x, y;
    curve_hilbert_d2xy(d, 8, &x, &y);

    uint32_t xs[4] = {0, 1, 2, 3}, ys[4] = {0, 0, 1, 1}, zs[4];
    curve_morton2_batch(xs, ys, zs, 4);
    curve_hilbert_xy2d_batch(xs, ys, 8, zs, 4);

    printf("%u %s\n", d, curve_compiled_isa());
    return 0;
}
```

```sh
gcc -std=c11 -O2 -I include main.c -o main
# x86: -mbmi2 -mavx2  (or -mavx512f)
# x86 fat binary: -DCURVE_RUNTIME_DISPATCH compiles extra AVX kernels
# aarch64: NEON batch kernels come along for the ride
```

That is the happy path. Everything else is either a wider key or the same
pattern with more axes.

## What it supports

| Family | Notes |
|---|---|
| Gray code | Binary-reflected, 64-bit |
| Morton / Z-order | 2D (16- and 32-bit/axis), 3D (21-bit/axis) |
| Hilbert 2D | Constant-time, bits 1..16; 64-bit variant bits 1..32 |
| Hilbert 3D | Morton + A4 lookup, bits 1..10 |
| Hilbert n-D | Skilling / Butz, `n_dims * bits ≤ 64` |
| Wide n-D keys | Checked Morton / Hilbert, up to 256-D × 32-bit, caller-owned words |
| Prefix / cells | Bit prefixes, parent/child, AABB bounds |
| Range queries | Inclusive AABB → merged `[lower, upper)` key spans |
| Prefix directory | Pointerless PATRICIA over sorted keys |
| Quantize / fused | Linear bins, fused float → Morton/Hilbert key |
| Moore | Cyclic Hilbert, bits 1..16 |
| Peano 2D | Geometric base-3, digits 1..10, invertible |
| Gray-coded Morton | Gray each axis, then interleave |
| Mixed-bit Morton | Anisotropic 1..32 bits per axis |
| Sierpiński | Ordering key, non-dense, xy2d only |
| Windowed | `double` coords through a rectangle onto a 16-bit grid |
| Gilbert 2D/3D | Visitor traversal of arbitrary rectangles/boxes |
| Batch | SoA; NEON / AVX2 / AVX-512 for Gray, Morton, Hilbert 2D |
| Runtime dispatch | Caller-owned `curve_dispatch`; opt-in `CURVE_RUNTIME_DISPATCH` |

## House rules

**Legacy maps clamp.** `bits` is the grid: side `2^bits`, index
`2^(dims*bits)`. Out-of-range `bits` on the original 32/64-bit APIs clamp.
They never assert. Coordinates on Hilbert / Moore / Sierpiński / n-D are
masked onto that grid. Extra high bits wrap. They do not pollute the index.

**Checked APIs return `curve_status`.** Wide keys, prefixes, range walks,
and quantizers reject bad arguments instead of clamping. Failed
encode/decode leaves the output buffer alone.

**`xy2d` goes out. `d2xy` comes back.** Decode pointers are `restrict`.
If they alias, that is your bug.

**Caller owns the scratch.** No hidden VLAs, no process globals, no
`malloc` inside the header. Hilbert wide-encode wants `dimensions`
`uint32_t`s. Range walks publish their scratch sizes. The vtable lives
in a `curve_dispatch` you pass in.

**SoA batches do not alias.** Unaligned pointers are fine. Define
`CURVE_NO_SIMD` to force scalar loops. BMI2 `PDEP`/`PEXT` still applies
to the scalar bit-interleave when you compiled with `-mbmi2`.

## Wide keys

`uint32_t` axes into a portable multi-word key. `dimensions` in
1..`CURVE_MAX_DIMS` (256), `bits_per_axis` in 1..32. Hilbert scratch may
alias `coords`. It must not alias the key.

```c
uint32_t coords[256], scratch[256];
uint64_t key[16];
uint8_t bytes[128];
curve_hilbertnd_encode_words(coords, 256, 4, key, 16, scratch, 256);
curve_key_words_to_bytes(key, 1024, bytes, 128);
```

`words[0]` holds the most significant 64 bits. Unused bits of the last
word sit in the low end and are zero. Bytes are the big-endian view of
those logical bits, host-endian independent. When `K = n_dims*bits ≤ 64`,
a wide Hilbert key matches the legacy integer after `words[0] >> (64-K)`.

Wide Morton uses the same plane order as n-D Hilbert (MSB level first,
axis 0 first within a level). That is a permutation of `curve_morton2` /
`curve_morton3`, not the same bit pattern. We documented it because the
other choice would have been a surprise later.

| `curve_status` | Meaning |
|---|---|
| `CURVE_OK` | success |
| `CURVE_EINVAL` | null, bad dims/bits/bounds, noncanonical tail, partial cell |
| `CURVE_ERANGE` | dims > max, no successor/parent, leaf has no child |
| `CURVE_ENOSPACE` | output or scratch too small |
| `CURVE_EOVERFLOW` | size arithmetic overflow |
| `CURVE_EUNSUPPORTED` | child index when `n_dims > 64`; box walk too wide |

## Prefixes, cells, ranges

Bit 0 of a wide key is its most significant logical bit. A prefix of
length `P` is those first `P` bits. The interval is `[masked, successor)`
when a successor exists, otherwise `[masked, +∞)`.

A **cell** is a prefix whose length is a multiple of `dimensions`. Partial
levels are `CURVE_EINVAL`. Child indices are `uint64_t` and need
`n_dims ≤ 64`. Bounds APIs work at any dimension count. Morton bounds
come off the interleaved prefix bits. Hilbert bounds decode a point in
the cell and snap to the dyadic hypercube.

```c
curve_cell cell, child;
uint64_t root = 0, child_key;
uint32_t mins[2], maxs[2], scratch[2];
curve_cell_init(&cell, 2, 8, 0);
curve_cell_child(&cell, &root, 0, &child, &child_key, 1);
curve_hilbert_cell_bounds(&child, &child_key, mins, maxs, 2, scratch, 2);
```

`curve_morton_box_intervals` / `curve_hilbert_box_intervals` walk cells
in key order and emit merged spans covering an inclusive integer AABB.
`n_dims` is 1..`CURVE_RANGE_MAX_DIMS` (8). Capacity 0 counts.
`CURVE_RANGE_EXACT` never includes a point outside the box.
`CURVE_RANGE_OVERAPPROX` never drops coverage. A pointerless PATRICIA
(`curve_prefix_dir_build`) sits on a sorted key array if you want one.

## Quantize

Linear per-axis maps from `float` / `double` onto `0 .. 2^bits - 1`.
No `<math.h>`. Compile without `-ffast-math` if you care about the same
answer on two machines.

| Mode | Use | Rule |
|---|---|---|
| `CURVE_QUANT_FLOOR` | geometric / GIS bins | equal-width intervals; `max` → last bin |
| `CURVE_QUANT_NEAREST` | embeddings | affine grid, `floor(s+0.5)`; ties toward +∞ |

NaN and −Inf go to 0. +Inf goes to the last bin. Finite underflow /
overflow saturate. A zero-width finite range maps everything to 0.
Fused `curve_mortonnd_encode_f64` / `curve_hilbertnd_encode_f64` quantize
into caller scratch, then pack a key. Bit-identical to the two-step.
Hilbert still materializes every axis. There is no free lunch in Skilling.

The older window helpers (`curve_z_window`, `curve_hilbert_window`) still
map a rectangle onto a 16-bit grid via `t * 0xFFFF`. They have not moved.

## Batch and dispatch

Public `curve_*_batch` functions pick the compile-time widest ISA.
`curve_compiled_isa()` (alias `curve_simd_isa()`) names that path.
32-bit Morton and Hilbert 2D use NEON / AVX2 / AVX-512. 64-bit Morton
and Hilbert 2D use AVX2 / AVX-512; 2-wide NEON lost to the scalar loop,
so it is not here. n-D Hilbert stays scalar. Skilling couples axes with
data-dependent swaps, and 256-D is that walk.

Runtime selection is a caller-owned `curve_dispatch`. No process globals.
`curve_dispatch_init(&d, allow)` fills pointers from compiled ∩ runtime ∩
`allow` (`allow == 0` means any). SCALAR is always available. x86 uses
CPUID + XGETBV. aarch64 NEON is baseline. Define `CURVE_RUNTIME_DISPATCH`
on x86 gcc/clang to compile extra AVX kernels with target attributes so a
baseline binary can still pick them up at run time.

```c
curve_dispatch d;
curve_dispatch_init(&d, 0);                 /* widest legal kernel */
d.morton2_batch(xs, ys, out, n);
curve_dispatch_init(&d, CURVE_ISA_SCALAR);  /* force scalar */
```

## A few curves, with the labels they earned

**Peano 2D** is the connected geometric curve on a `3^d × 3^d` grid
(digits 1..10, index `0 .. 9^d - 1`). Consecutive cells share an edge.
This is not ternary Morton. If you wanted digit interleave, you already
have Morton.

**Gray-coded Morton** Gray-encodes each axis after masking, then
interleaves. Inverse compact then Gray-decode.

**Mixed-bit Morton** packs anisotropic axes: plane `p` from
`max(bits)-1` down to 0 emits axis `i` iff `bits[i] > p`. Uniform widths
match `curve_mortonnd_encode_words`. The cells are rectangles. They are
not Hilbert orthants, and calling them that would be a lie.

**Moore** is Hilbert with the last cell glued to the first. Order 1
equals order-1 Hilbert.

**Sierpiński** is an ordering key, not a dense index. Keys live in
`[0, 2^(2*bits+1))` with half the space unused. Sorting by key yields
the curve (Bartholdi & Platzman 1989). No `d2xy`.

**Hilbert n-D** takes a coordinate array. `n_dims < 1` is a no-op. n=2
matches `curve_hilbert_xy2d`. n=3 does not match `curve_hilbert3_xyz2d`.
Both are valid Hilbert curves. They just do not share an origin story.

**Gilbert** fills an arbitrary rectangle via a visitor. 2D paths are
edge-adjacent. 3D paths are face-adjacent when every side is even. Odd
sides may step diagonally between sub-blocks, as in the reference. That
is a known upstream limitation, not a transcription bug.

## Build and tests

You need a C99 compiler. CMake 3.14 or newer is optional. None of this
is exotic. If it fails to build, the problem is probably local.

```sh
gcc -std=c11 -pedantic -Wall -Wextra -O2 -I include tests/test_curves.c -lm -o test_curves
./test_curves

cmake -S . -B build -DCURVES_BUILD_BENCHMARKS=ON
cmake --build build
ctest --test-dir build --output-on-failure
./build/benchmarks/bench_main
```

CI runs gcc/clang × ASan+UBSan, C99 pedantic, `CURVE_NO_SIMD`, and
`CURVE_RUNTIME_DISPATCH`. The suite also compiles as C++11. Hot-path
SIMD output is required to match scalar byte for byte.

| CMake option | Default | Effect |
|---|---|---|
| `CURVES_BUILD_TESTS` | ON | C suite + C++11 smoke test |
| `CURVES_BUILD_EXAMPLES` | ON | The programs in `examples/` |
| `CURVES_BUILD_BENCHMARKS` | OFF | `bench_main` (`-O3 -march=native`) |
| `CURVES_ENABLE_SANITIZERS` | OFF | ASan + UBSan, no recover |

```cmake
find_package(curves REQUIRED)
target_link_libraries(your_target PRIVATE curves::curves)
```

or `add_subdirectory` and link `curves::curves` the same way. The
`curves::curves` target is INTERFACE: include path, C11, nothing to
link against but libm for the tests.

## Examples

Each file in `examples/` is one program you can read top to bottom.
They compile. They are starting points.

| Program | What you get |
|---|---|
| `render_pgm` | Morton / Hilbert / Moore / Sierpiński / Gilbert as `.pgm` |
| `locality_demo` | Row-major vs Morton vs Hilbert vs Moore |
| `wide_key_demo` | 256-D × 4-bit Morton and Hilbert keys |
| `prefix_cell_demo` | Prefix range + Hilbert cell AABB |
| `range_query_demo` | AABB → Morton/Hilbert key intervals |
| `quantize_demo` | float → 4-bit bins → 256-D keys |
| `batch_demo` | Hilbert 2D SIMD batch vs scalar |
| `research_demo` | Peano, Gray-Morton, mixed-bit, dispatch |

```sh
./build/examples/render_pgm /tmp/out
./build/examples/research_demo
```

## Performance

Numbers from an ARM Cortex-A78AE, Clang 22, `-O3 -march=native`. If a
SIMD kernel loses to the scalar loop, it is not in the header. See
`benchmarks/results/`.

| Operation | ns/op | Mops/s |
|---|---|---|
| `gray_encode` | 3.4 | 292 |
| `morton2` encode | 4.3 | 233 |
| `morton2` batch | 3.0 | 337 |
| `hilbert2` xy2d (16-bit) | 18 | 54 |
| `hilbert2` xy2d batch (NEON) | 11 | 95 |
| `mortonnd` 256-D × 4-bit | 3600 | 0.28 |
| `hilbertnd` 256-D × 4-bit | 8000 | 0.13 |

## Provenance

- 2D Hilbert and windowed maps: [rawrunprotected/hilbert_curves](https://github.com/rawrunprotected/hilbert_curves) and [tidwall/curve](https://github.com/tidwall/curve) (public domain)
- 3D Hilbert A4 tables: same source
- Gilbert: [jakubcerveny/gilbert](https://github.com/jakubcerveny/gilbert) (BSD-2-Clause) via skookum/spatial
- n-D Hilbert: A. R. Butz (1971); Doug Moore's Rice implementation; formulation after J. Skilling (2004)
- Moore: E. H. Moore (1900); quadrant structure from its L-system
- Sierpiński ordering: Bartholdi & Platzman (1989)
- Peano 2D: G. Peano (1890); this header uses the connected 9-way generator, not ternary-digit interleave

The comments at the top of `include/curves.h` keep the attribution next
to the code. Leave them there.

## License

MIT. See [LICENSE](LICENSE). Keep the notice with redistributed copies.
The code is more useful when nobody has to guess what they are allowed
to do with it. If you improve it, send the patch. If you do not, that
is also fine.
