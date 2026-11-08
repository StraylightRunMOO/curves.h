# curves.h 1.0.0 baseline — aarch64

Recorded before wide-key work on `feature/wide-spatial-keys`.

| Field | Value |
|---|---|
| Date | 2026-11-07 |
| Host | ARM Cortex-A78AE, 6 cores, little-endian aarch64 |
| ISA (`curve_simd_isa()`) | `neon` |
| BMI2 | portable fallback (not an x86 host) |
| CMake compiler | Clang 22.1.3 |
| Flags | `-O3 -march=native` (benchmarks), `-O2` (tests) |
| Tests | 78 passed, 0 failed (`test_curves` at commit `22e4840`); extra adjacency/reference tests added in PR-A |

Human-readable `bench_main` (Release, Clang 22, `-O3 -march=native`):

```
gray_encode                    269.40 Mops/s  (    3.71 ns/op)
gray_decode                    286.37 Mops/s  (    3.49 ns/op)
morton2 encode                 233.75 Mops/s  (    4.28 ns/op)
morton2 decode                 219.99 Mops/s  (    4.55 ns/op)
morton2_64 encode               93.77 Mops/s  (   10.66 ns/op)
morton2_64 decode              198.89 Mops/s  (    5.03 ns/op)
morton3 encode                 243.03 Mops/s  (    4.11 ns/op)
morton3 decode                 137.15 Mops/s  (    7.29 ns/op)
hilbert2 xy2d (16-bit)          54.25 Mops/s  (   18.43 ns/op)
hilbert2 d2xy (16-bit)         100.76 Mops/s  (    9.93 ns/op)
hilbert2_64 xy2d (32-bit)       33.71 Mops/s  (   29.67 ns/op)
hilbert2_64 d2xy (32-bit)       84.99 Mops/s  (   11.77 ns/op)
hilbert3 xyz2d (10-bit)         37.96 Mops/s  (   26.34 ns/op)
hilbert3 d2xyz (10-bit)         34.96 Mops/s  (   28.60 ns/op)
hilbertnd xy2d (n=3,b=10)       11.50 Mops/s  (   86.92 ns/op)
hilbertnd d2xy (n=3,b=10)       16.53 Mops/s  (   60.48 ns/op)
moore xy2d (16-bit)             47.85 Mops/s  (   20.90 ns/op)
moore d2xy (16-bit)             84.00 Mops/s  (   11.91 ns/op)
sierpinski xy2d (16-bit)        12.73 Mops/s  (   78.54 ns/op)
morton2 batch                  401.76 Mops/s  (    2.49 ns/op)
morton2 decode batch           401.59 Mops/s  (    2.49 ns/op)
morton3 batch                  130.82 Mops/s  (    7.64 ns/op)
hilbert2 xy2d batch            106.61 Mops/s  (    9.38 ns/op)
gilbert2 1920x1080             127.31 Mcells/s
gilbert3 64x64x64               88.86 Mcells/s
```

These are single-shot wall times on a loaded workstation, not a controlled performance lab. Use them as a regression fence (~10%) for existing scalar APIs, not as published throughput claims.

CSV form: `./build/benchmarks/bench_main --csv`.

## v1.1.0 (same host, after wide keys)

Existing scalar APIs are within a few percent of the 1.0.0 fence (no >10%
regression). New wide-key rows from `v1.1-aarch64.csv` (Clang 22, `-O3
-march=native`, inputs vary per iteration so the encoder is not CSE'd):

| Operation | ns/key | notes |
|---|---|---|
| mortonnd enc n=2 b=16 | ~21 | scalar bit loop; slower than SWAR `curve_morton2` |
| hilbertnd enc n=2 b=16 | ~101 | vs legacy `hilbertnd xy2d (n=3,b=10)` ~87 ns |
| mortonnd enc n=256 b=4 | ~2990 | 1024-bit key, 16 words, 1 KB coords |
| hilbertnd enc n=256 b=4 | ~6460 | Skilling + pack; scratch = 1 KB |
| mortonnd enc n=256 b=8 | ~5980 | 2048-bit key |
| hilbertnd enc n=256 b=8 | ~13620 | |

These are single-shot wall times, not a controlled lab. Scratch for 256-D
Hilbert is 256 × 4 = 1024 bytes, caller-owned.
