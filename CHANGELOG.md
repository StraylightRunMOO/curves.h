# Changelog

All notable changes to this project will be documented in this file.

## [1.1.0] - 2026-11-07

Spatial keys, hierarchy, SIMD Hilbert, and a few more curves. The 1.0
32/64-bit maps stay bit-identical.

### Added

- Checked wide n-D Morton / Hilbert keys (`curve_status`, caller-owned
  words, up to 256-D × 32-bit). `words[0]` is the MSBs; unused tail bits
  are left-aligned. Bytes are big-endian of those logical bits.
- Prefix ops, spatial cells (parent/child/sibling, Morton/Hilbert AABB
  bounds), AABB → merged key-interval decomposition, and a pointerless
  PATRICIA directory over sorted keys.
- Linear quantizers (`CURVE_QUANT_FLOOR` / `CURVE_QUANT_NEAREST`) and
  fused float → Morton/Hilbert encode. No `<math.h>`.
- SIMD Hilbert 2D batches (NEON / AVX2 / AVX-512). `curve_compiled_isa()`
  names the compile-time path. 64-bit Hilbert batch is AVX2/AVX-512 only.
- Caller-owned `curve_dispatch` vtable. Opt-in `CURVE_RUNTIME_DISPATCH`
  compiles extra x86 AVX kernels with target attributes. No process globals.
- Gray-coded Morton 2D/3D/n-D, geometric Peano 2D (digits 1..10,
  invertible), mixed-bit (anisotropic) Morton wide keys.

### Notes

Wide Morton packing matches n-D Hilbert (MSB level first, axis 0 first
within a level). That is a permutation of `curve_morton2` / `curve_morton3`.
Checked APIs reject invalid arguments. Legacy functions still clamp.

Not in 1.1: GPU kernels, learned embeddings, prefix-compressed persistence,
mixed-bit Hilbert, Peano 3D, n-D Hilbert SIMD.

## [1.0.0] - 2026-08-23

Initial release. Header-only C99/C11 Gray, Morton 2D/3D, Hilbert 2D/3D/n-D,
Moore, Sierpiński, windowed maps, Gilbert 2D/3D, and SIMD batch APIs for
Gray and Morton.
