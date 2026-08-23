/* ============================================================================
 * curves.h — Space-filling curves, Gray codes, and spatial orderings
 *
 * Single-file, header-only C99/C11 library of locality-preserving maps:
 *
 *   • Binary-reflected Gray encode / decode
 *   • Morton (Z-order) 2D / 3D, 32- and 64-bit  +  windowed Z-curve
 *   • Hilbert 2D (constant-time, bits 1..16) and 64-bit variant (bits 1..32)
 *     + windowed Hilbert
 *   • Hilbert 3D (Morton + A4 lookup, bits 1..10)
 *   • Hilbert n-D (Skilling formulation of the Butz algorithm, n*bits <= 64)
 *   • Moore curve (cyclic Hilbert variant, bits 1..16)
 *   • Sierpiński curve ordering key (Bartholdi & Platzman; non-dense, xy2d only)
 *   • Gilbert 2D / 3D visit (arbitrary rectangles, not just 2^n)
 *   • Batch APIs (restrict + NEON / AVX2 / AVX-512) for Gray, Morton, Hilbert
 *
 * 2D Hilbert / windowed maps follow:
 *   rawrunprotected/hilbert_curves  (public domain)
 *   tidwall/curve                   (public domain)
 * 3D Hilbert lookup tables: same source.
 * Gilbert: jakubcerveny/gilbert via skookum/spatial.
 * n-D Hilbert: A. R. Butz (1971); Doug Moore's Rice U. implementation;
 *   formulation after J. Skilling, "Programming the Hilbert curve" (2004).
 * Moore: E. H. Moore (1900); quadrant structure derived from its L-system.
 * Sierpiński ordering: Bartholdi & Platzman, spacefilling-curve TSP (1989).
 *
 * Compile: -std=c11 -pedantic -Wall -Wextra
 * Optional: -mbmi2 (x86), -mavx2 / -mavx512f, ARM NEON is used automatically
 *   on aarch64.  Define CURVE_NO_SIMD to force scalar batch loops.
 *
 * Version: 1.0.0
 * License: MIT
 * ============================================================================ */

#ifndef CURVES_H
#define CURVES_H

#include <stdint.h>
#include <stddef.h>

#define CURVES_VERSION_MAJOR 1
#define CURVES_VERSION_MINOR 0
#define CURVES_VERSION_PATCH 0
#define CURVES_VERSION       "1.0.0"

/* -------------------------------------------------------------------------- */
/* Compiler / ISA helpers                                                     */
/* -------------------------------------------------------------------------- */

#if defined(__AVX512F__) || defined(__AVX2__) || defined(__BMI2__)
#  include <immintrin.h>
#endif
#if defined(__BMI2__)
#  define CURVE_HAVE_BMI2 1
#endif
#if !defined(CURVE_NO_SIMD)
#  if defined(__AVX512F__)
#    define CURVE_HAVE_AVX512 1
#  endif
#  if defined(__AVX2__)
#    define CURVE_HAVE_AVX2 1
#  endif
#  if defined(__ARM_NEON) || defined(__ARM_NEON__) || defined(__aarch64__)
#    include <arm_neon.h>
#    define CURVE_HAVE_NEON 1
#  endif
#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifndef CURVE_INLINE
#  define CURVE_INLINE static inline
#endif

/* C99 `restrict`; GNU/MSVC `__restrict` in C++. */
#ifndef CURVE_RESTRICT
#  ifdef __cplusplus
#    if defined(__GNUC__) || defined(__clang__) || defined(_MSC_VER)
#      define CURVE_RESTRICT __restrict
#    else
#      define CURVE_RESTRICT
#    endif
#  else
#    define CURVE_RESTRICT restrict
#  endif
#endif

#ifndef CURVE_ALIGNAS
#  ifdef __cplusplus
#    define CURVE_ALIGNAS(n) alignas(n)
#  elif defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L
#    define CURVE_ALIGNAS(n) _Alignas(n)
#  elif defined(_MSC_VER)
#    define CURVE_ALIGNAS(n) __declspec(align(n))
#  elif defined(__GNUC__)
#    define CURVE_ALIGNAS(n) __attribute__((aligned(n)))
#  else
#    define CURVE_ALIGNAS(n)
#  endif
#endif

#ifdef __cplusplus
#  define CURVE_STATIC_ASSERT(cond, msg) static_assert(cond, msg)
#elif defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L
#  define CURVE_STATIC_ASSERT(cond, msg) _Static_assert(cond, msg)
#else
#  define CURVE_STATIC_ASSERT(cond, msg)
#endif

#if defined(__GNUC__) || defined(__clang__)
#  define CURVE_LIKELY(x)    __builtin_expect(!!(x), 1)
#  define CURVE_UNLIKELY(x)  __builtin_expect(!!(x), 0)
#  define CURVE_PREFETCH(p)  __builtin_prefetch((const void *)(p), 0, 3)
#  define CURVE_ASSUME(x)    do { if (!(x)) __builtin_unreachable(); } while (0)
#else
#  define CURVE_LIKELY(x)    (x)
#  define CURVE_UNLIKELY(x)  (x)
#  define CURVE_PREFETCH(p)  ((void)0)
#  define CURVE_ASSUME(x)    ((void)0)
#endif

/* Low `bits` bits set; well-defined for bits == 32 (unlike `1u << bits`). */
CURVE_INLINE uint32_t curve__lo_mask(int bits)
{
    return (bits >= 32) ? 0xFFFFFFFFu : ((1u << bits) - 1u);
}

/* -------------------------------------------------------------------------- */
/* Gray codes (binary reflected)                                              */
/* -------------------------------------------------------------------------- */

CURVE_INLINE uint64_t curve_gray_encode(uint64_t x)
{
    return x ^ (x >> 1);
}

CURVE_INLINE uint64_t curve_gray_decode(uint64_t g)
{
    g ^= g >> 32;
    g ^= g >> 16;
    g ^= g >> 8;
    g ^= g >> 4;
    g ^= g >> 2;
    g ^= g >> 1;
    return g;
}

/* -------------------------------------------------------------------------- */
/* Morton / Z-order                                                           */
/* -------------------------------------------------------------------------- */

CURVE_INLINE uint32_t curve__part1by1(uint32_t x)
{
#if defined(CURVE_HAVE_BMI2)
    return _pdep_u32(x, 0x55555555u);
#else
    x &= 0x0000ffffu;
    x = (x | (x << 8)) & 0x00ff00ffu;
    x = (x | (x << 4)) & 0x0f0f0f0fu;
    x = (x | (x << 2)) & 0x33333333u;
    x = (x | (x << 1)) & 0x55555555u;
    return x;
#endif
}

CURVE_INLINE uint32_t curve__compact1by1(uint32_t x)
{
#if defined(CURVE_HAVE_BMI2)
    return _pext_u32(x, 0x55555555u);
#else
    x &= 0x55555555u;
    x = (x | (x >> 1)) & 0x33333333u;
    x = (x | (x >> 2)) & 0x0f0f0f0fu;
    x = (x | (x >> 4)) & 0x00ff00ffu;
    x = (x | (x >> 8)) & 0x0000ffffu;
    return x;
#endif
}

CURVE_INLINE uint32_t curve_morton2(uint32_t x, uint32_t y)
{
    return curve__part1by1(x) | (curve__part1by1(y) << 1);
}

CURVE_INLINE void curve_morton2_decode(uint32_t d, uint32_t *CURVE_RESTRICT x,
                                       uint32_t *CURVE_RESTRICT y)
{
    *x = curve__compact1by1(d);
    *y = curve__compact1by1(d >> 1);
}

/* 64-bit Morton: 32 bits per axis interleaved into uint64. */
CURVE_INLINE uint64_t curve__part1by1_64(uint64_t x)
{
#if defined(CURVE_HAVE_BMI2)
    return _pdep_u64(x, 0x5555555555555555ull);
#else
    x &= 0x00000000ffffffffull;
    x = (x | (x << 16)) & 0x0000ffff0000ffffull;
    x = (x | (x <<  8)) & 0x00ff00ff00ff00ffull;
    x = (x | (x <<  4)) & 0x0f0f0f0f0f0f0f0full;
    x = (x | (x <<  2)) & 0x3333333333333333ull;
    x = (x | (x <<  1)) & 0x5555555555555555ull;
    return x;
#endif
}

CURVE_INLINE uint64_t curve__compact1by1_64(uint64_t x)
{
#if defined(CURVE_HAVE_BMI2)
    return _pext_u64(x, 0x5555555555555555ull);
#else
    x &= 0x5555555555555555ull;
    x = (x | (x >>  1)) & 0x3333333333333333ull;
    x = (x | (x >>  2)) & 0x0f0f0f0f0f0f0f0full;
    x = (x | (x >>  4)) & 0x00ff00ff00ff00ffull;
    x = (x | (x >>  8)) & 0x0000ffff0000ffffull;
    x = (x | (x >> 16)) & 0x00000000ffffffffull;
    return x;
#endif
}

CURVE_INLINE uint64_t curve_morton2_64(uint32_t x, uint32_t y)
{
    return curve__part1by1_64(x) | (curve__part1by1_64(y) << 1);
}

CURVE_INLINE void curve_morton2_64_decode(uint64_t d, uint32_t *CURVE_RESTRICT x,
                                          uint32_t *CURVE_RESTRICT y)
{
    *x = (uint32_t)curve__compact1by1_64(d);
    *y = (uint32_t)curve__compact1by1_64(d >> 1);
}

CURVE_INLINE uint64_t curve__part1by2(uint64_t x)
{
#if defined(CURVE_HAVE_BMI2)
    return _pdep_u64(x, 0x1249249249249249ull);
#else
    x &= 0x1fffffull;
    x = (x | (x << 32)) & 0x1f00000000ffffull;
    x = (x | (x << 16)) & 0x1f0000ff0000ffull;
    x = (x | (x <<  8)) & 0x100f00f00f00f00full;
    x = (x | (x <<  4)) & 0x10c30c30c30c30c3ull;
    x = (x | (x <<  2)) & 0x1249249249249249ull;
    return x;
#endif
}

CURVE_INLINE uint64_t curve__compact1by2(uint64_t x)
{
#if defined(CURVE_HAVE_BMI2)
    return _pext_u64(x, 0x1249249249249249ull);
#else
    x &= 0x1249249249249249ull;
    x = (x | (x >>  2)) & 0x10c30c30c30c30c3ull;
    x = (x | (x >>  4)) & 0x100f00f00f00f00full;
    x = (x | (x >>  8)) & 0x1f0000ff0000ffull;
    x = (x | (x >> 16)) & 0x1f00000000ffffull;
    x = (x | (x >> 32)) & 0x1fffffull;
    return x;
#endif
}

CURVE_INLINE uint64_t curve_morton3(uint32_t x, uint32_t y, uint32_t z)
{
    return curve__part1by2(x) | (curve__part1by2(y) << 1) | (curve__part1by2(z) << 2);
}

CURVE_INLINE void curve_morton3_decode(uint64_t d, uint32_t *CURVE_RESTRICT x,
                                       uint32_t *CURVE_RESTRICT y, uint32_t *CURVE_RESTRICT z)
{
    *x = (uint32_t)curve__compact1by2(d);
    *y = (uint32_t)curve__compact1by2(d >> 1);
    *z = (uint32_t)curve__compact1by2(d >> 2);
}

/* 10-bit-per-axis Morton packed into uint32 (used by 3D Hilbert). */
CURVE_INLINE uint32_t curve__part1by2_u32(uint32_t x)
{
#if defined(CURVE_HAVE_BMI2)
    return _pdep_u32(x, 0x09249249u);
#else
    x &= 0x3ffu;
    x = (x | (x << 16)) & 0x030000FFu;
    x = (x | (x <<  8)) & 0x0300F00Fu;
    x = (x | (x <<  4)) & 0x030C30C3u;
    x = (x | (x <<  2)) & 0x09249249u;
    return x;
#endif
}

CURVE_INLINE uint32_t curve__compact1by2_u32(uint32_t x)
{
#if defined(CURVE_HAVE_BMI2)
    return _pext_u32(x, 0x09249249u);
#else
    x &= 0x09249249u;
    x = (x | (x >>  2)) & 0x030C30C3u;
    x = (x | (x >>  4)) & 0x0300F00Fu;
    x = (x | (x >>  8)) & 0x030000FFu;
    x = (x | (x >> 16)) & 0x000003FFu;
    return x;
#endif
}

/* -------------------------------------------------------------------------- */
/* Window helpers (tidwall/curve)                                             */
/* -------------------------------------------------------------------------- */

CURVE_INLINE uint16_t curve__unit_u16(double t)
{
    /* NaN and -Inf fail `t >= 0`; +Inf is clamped by the second test. */
    if (!(t >= 0.0)) t = 0.0;
    if (t > 1.0) t = 1.0;
    return (uint16_t)(t * 0xFFFF);
}

CURVE_INLINE void curve__window_xy(double x, double y,
                                   const double *CURVE_RESTRICT window,
                                   uint16_t *CURVE_RESTRICT xi,
                                   uint16_t *CURVE_RESTRICT yi)
{
    double dx = window[2] - window[0];
    double dy = window[3] - window[1];
    *xi = curve__unit_u16(dx != 0.0 ? (x - window[0]) / dx : 0.0);
    *yi = curve__unit_u16(dy != 0.0 ? (y - window[1]) / dy : 0.0);
}

/* Windowed Z-order.  window = {xmin, ymin, xmax, ymax}. */
CURVE_INLINE uint32_t curve_z_window(double x, double y,
                                     const double *CURVE_RESTRICT window)
{
    uint16_t xint, yint;
    curve__window_xy(x, y, window, &xint, &yint);
    return curve_morton2((uint32_t)xint, (uint32_t)yint);
}

/* -------------------------------------------------------------------------- */
/* 2D Hilbert — constant-time (rawrunprotected / tidwall), bits in 1..16      */
/* -------------------------------------------------------------------------- */

CURVE_INLINE int curve__clamp_bits2(int bits)
{
    if (bits < 1)  return 1;
    if (bits > 16) return 16;
    return bits;
}

CURVE_INLINE uint32_t curve__prefixscan(uint32_t x)
{
    x = (x >> 8) ^ x;
    x = (x >> 4) ^ x;
    x = (x >> 2) ^ x;
    x = (x >> 1) ^ x;
    return x;
}

CURVE_INLINE uint32_t curve_hilbert_xy2d(uint32_t x, uint32_t y, int bits)
{
    bits = curve__clamp_bits2(bits);
    x &= curve__lo_mask(bits);
    y &= curve__lo_mask(bits);
    x <<= (16 - bits);
    y <<= (16 - bits);

    uint32_t a = x ^ y;
    uint32_t b = 0xFFFFu ^ a;
    uint32_t c = 0xFFFFu ^ (x | y);
    uint32_t d = x & (y ^ 0xFFFFu);

    uint32_t A = a | (b >> 1);
    uint32_t B = (a >> 1) ^ a;
    uint32_t C = ((c >> 1) ^ (b & (d >> 1))) ^ c;
    uint32_t D = ((a & (c >> 1)) ^ (d >> 1)) ^ d;

    a = A; b = B; c = C; d = D;
    A = ((a & (a >> 2)) ^ (b & (b >> 2)));
    B = ((a & (b >> 2)) ^ (b & ((a ^ b) >> 2)));
    C ^= ((a & (c >> 2)) ^ (b & (d >> 2)));
    D ^= ((b & (c >> 2)) ^ ((a ^ b) & (d >> 2)));

    a = A; b = B; c = C; d = D;
    A = ((a & (a >> 4)) ^ (b & (b >> 4)));
    B = ((a & (b >> 4)) ^ (b & ((a ^ b) >> 4)));
    C ^= ((a & (c >> 4)) ^ (b & (d >> 4)));
    D ^= ((b & (c >> 4)) ^ ((a ^ b) & (d >> 4)));

    a = A; b = B; c = C; d = D;
    C ^= ((a & (c >> 8)) ^ (b & (d >> 8)));
    D ^= ((b & (c >> 8)) ^ ((a ^ b) & (d >> 8)));

    a = C ^ (C >> 1);
    b = D ^ (D >> 1);

    uint32_t i0 = x ^ y;
    uint32_t i1 = b | (0xFFFFu ^ (i0 | a));

    return (curve__part1by1(i1) << 1 | curve__part1by1(i0)) >> (32 - 2 * bits);
}

CURVE_INLINE void curve_hilbert_d2xy(uint32_t d, int bits,
                                     uint32_t *CURVE_RESTRICT x,
                                     uint32_t *CURVE_RESTRICT y)
{
    bits = curve__clamp_bits2(bits);
    uint32_t i = d << (32 - 2 * bits);

    uint32_t i0 = curve__compact1by1(i);
    uint32_t i1 = curve__compact1by1(i >> 1);

    uint32_t t0 = (i0 | i1) ^ 0xFFFFu;
    uint32_t t1 = i0 & i1;

    uint32_t a = (((i0 ^ 0xFFFFu) & curve__prefixscan(t1)) |
                  (i0 & curve__prefixscan(t0)));

    *x = (a ^ i1) >> (16 - bits);
    *y = (a ^ i0 ^ i1) >> (16 - bits);
}

/* Windowed Hilbert.  window = {xmin, ymin, xmax, ymax}.  Always 16-bit. */
CURVE_INLINE uint32_t curve_hilbert_window(double x, double y,
                                          const double *CURVE_RESTRICT window)
{
    uint16_t xint, yint;
    curve__window_xy(x, y, window, &xint, &yint);
    return curve_hilbert_xy2d((uint32_t)xint, (uint32_t)yint, 16);
}

/* -------------------------------------------------------------------------- */
/* 2D Hilbert, 64-bit index — same algorithm widened, bits in 1..32           */
/* -------------------------------------------------------------------------- */

CURVE_INLINE int curve__clamp_bits2_64(int bits)
{
    if (bits < 1)  return 1;
    if (bits > 32) return 32;
    return bits;
}

CURVE_INLINE uint64_t curve__prefixscan64(uint64_t x)
{
    x = (x >> 16) ^ x;
    x = (x >>  8) ^ x;
    x = (x >>  4) ^ x;
    x = (x >>  2) ^ x;
    x = (x >>  1) ^ x;
    return x;
}

CURVE_INLINE uint64_t curve_hilbert2_64_xy2d(uint32_t x, uint32_t y, int bits)
{
    const uint64_t M = 0xFFFFFFFFull;
    uint64_t a, b, c, d, A, B, C, D, i0, i1;

    bits = curve__clamp_bits2_64(bits);
    x &= curve__lo_mask(bits);
    y &= curve__lo_mask(bits);
    uint64_t X = (uint64_t)x << (32 - bits);
    uint64_t Y = (uint64_t)y << (32 - bits);

    a = X ^ Y;
    b = M ^ a;
    c = M ^ (X | Y);
    d = X & (Y ^ M);

    A = a | (b >> 1);
    B = (a >> 1) ^ a;
    C = ((c >> 1) ^ (b & (d >> 1))) ^ c;
    D = ((a & (c >> 1)) ^ (d >> 1)) ^ d;

    a = A; b = B; c = C; d = D;
    A = ((a & (a >> 2)) ^ (b & (b >> 2)));
    B = ((a & (b >> 2)) ^ (b & ((a ^ b) >> 2)));
    C ^= ((a & (c >> 2)) ^ (b & (d >> 2)));
    D ^= ((b & (c >> 2)) ^ ((a ^ b) & (d >> 2)));

    a = A; b = B; c = C; d = D;
    A = ((a & (a >> 4)) ^ (b & (b >> 4)));
    B = ((a & (b >> 4)) ^ (b & ((a ^ b) >> 4)));
    C ^= ((a & (c >> 4)) ^ (b & (d >> 4)));
    D ^= ((b & (c >> 4)) ^ ((a ^ b) & (d >> 4)));

    a = A; b = B; c = C; d = D;
    A = ((a & (a >> 8)) ^ (b & (b >> 8)));
    B = ((a & (b >> 8)) ^ (b & ((a ^ b) >> 8)));
    C ^= ((a & (c >> 8)) ^ (b & (d >> 8)));
    D ^= ((b & (c >> 8)) ^ ((a ^ b) & (d >> 8)));

    a = A; b = B; c = C; d = D;
    C ^= ((a & (c >> 16)) ^ (b & (d >> 16)));
    D ^= ((b & (c >> 16)) ^ ((a ^ b) & (d >> 16)));

    a = C ^ (C >> 1);
    b = D ^ (D >> 1);

    i0 = X ^ Y;
    i1 = b | (M ^ (i0 | a));

    return (curve__part1by1_64(i1) << 1 | curve__part1by1_64(i0)) >> (64 - 2 * bits);
}

CURVE_INLINE void curve_hilbert2_64_d2xy(uint64_t d, int bits,
                                         uint32_t *CURVE_RESTRICT x,
                                         uint32_t *CURVE_RESTRICT y)
{
    const uint64_t M = 0xFFFFFFFFull;
    uint64_t i, i0, i1, t0, t1, a;

    bits = curve__clamp_bits2_64(bits);
    i = d << (64 - 2 * bits);

    i0 = curve__compact1by1_64(i);
    i1 = curve__compact1by1_64(i >> 1);

    t0 = (i0 | i1) ^ M;
    t1 = i0 & i1;

    a = (((i0 ^ M) & curve__prefixscan64(t1)) |
         (i0 & curve__prefixscan64(t0)));

    *x = (uint32_t)((a ^ i1) >> (32 - bits));
    *y = (uint32_t)((a ^ i0 ^ i1) >> (32 - bits));
}

/* -------------------------------------------------------------------------- */
/* Moore curve — cyclic Hilbert variant (E. H. Moore, 1900), bits in 1..16    */
/*                                                                            */
/* Order n >= 2: four Hilbert(n-1) quadrants rotated so the path closes into  */
/* a loop (derived from the Moore L-system).  Quadrant q = top 2 bits of d:   */
/*   q=0: bottom-left, rot90     q=1: top-left, rot90                         */
/*   q=2: top-right,  rot270     q=3: bottom-right, rot270                    */
/* Order 1 is identical to order-1 Hilbert.                                   */
/* -------------------------------------------------------------------------- */

CURVE_INLINE uint32_t curve_moore_xy2d(uint32_t x, uint32_t y, int bits)
{
    uint32_t h, hx, hy, q;

    bits = curve__clamp_bits2(bits);
    x &= curve__lo_mask(bits);
    y &= curve__lo_mask(bits);
    if (bits == 1) return curve_hilbert_xy2d(x, y, 1);

    h = 1u << (bits - 1);
    if (x < h && y < h)       { q = 0; hx = y;             hy = h - 1 - x; }
    else if (x < h)           { q = 1; hx = y - h;         hy = h - 1 - x; }
    else if (y >= h)          { q = 2; hx = h - 1 - (y - h); hy = x - h; }
    else                      { q = 3; hx = h - 1 - y;     hy = x - h; }

    return (q << (2 * bits - 2)) | curve_hilbert_xy2d(hx, hy, bits - 1);
}

CURVE_INLINE void curve_moore_d2xy(uint32_t d, int bits,
                                   uint32_t *CURVE_RESTRICT x,
                                   uint32_t *CURVE_RESTRICT y)
{
    uint32_t h, q, r, hx, hy;

    bits = curve__clamp_bits2(bits);
    if (bits == 1) { curve_hilbert_d2xy(d, 1, x, y); return; }

    h = 1u << (bits - 1);
    q = d >> (2 * bits - 2);
    r = d & ((1u << (2 * bits - 2)) - 1u);
    curve_hilbert_d2xy(r, bits - 1, &hx, &hy);

    switch (q) {
    case 0:  *x = h - 1 - hy; *y = hx;             break; /* rot90,  BL */
    case 1:  *x = h - 1 - hy; *y = h + hx;         break; /* rot90,  TL */
    case 2:  *x = h + hy;     *y = 2 * h - 1 - hx; break; /* rot270, TR */
    default: *x = h + hy;     *y = h - 1 - hx;     break; /* rot270, BR */
    }
}

/* -------------------------------------------------------------------------- */
/* 3D Hilbert — Morton + A4 lookup (rawrunprotected), bits in 1..10           */
/* -------------------------------------------------------------------------- */

CURVE_ALIGNAS(64) static const uint8_t curve__morton_to_hilbert3[] = {
    48, 33, 27, 34, 47, 78, 28, 77,
    66, 29, 51, 52, 65, 30, 72, 63,
    76, 95, 75, 24, 53, 54, 82, 81,
    18,  3, 17, 80, 61,  4, 62, 15,
     0, 59, 71, 60, 49, 50, 86, 85,
    84, 83,  5, 90, 79, 56,  6, 89,
    32, 23,  1, 94, 11, 12,  2, 93,
    42, 41, 13, 14, 35, 88, 36, 31,
    92, 37, 87, 38, 91, 74,  8, 73,
    46, 45,  9, 10,  7, 20, 64, 19,
    70, 25, 39, 16, 69, 26, 44, 43,
    22, 55, 21, 68, 57, 40, 58, 67,
};

CURVE_ALIGNAS(64) static const uint8_t curve__hilbert_to_morton3[] = {
    48, 33, 35, 26, 30, 79, 77, 44,
    78, 68, 64, 50, 51, 25, 29, 63,
    27, 87, 86, 74, 72, 52, 53, 89,
    83, 18, 16,  1,  5, 60, 62, 15,
     0, 52, 53, 57, 59, 87, 86, 66,
    61, 95, 91, 81, 80,  2,  6, 76,
    32,  2,  6, 12, 13, 95, 91, 17,
    93, 41, 40, 36, 38, 10, 11, 31,
    14, 79, 77, 92, 88, 33, 35, 82,
    70, 10, 11, 23, 21, 41, 40,  4,
    19, 25, 29, 47, 46, 68, 64, 34,
    45, 60, 62, 71, 67, 18, 16, 49,
};

CURVE_INLINE int curve__clamp_bits3(int bits)
{
    if (bits < 1)  return 1;
    if (bits > 10) return 10;
    return bits;
}

CURVE_STATIC_ASSERT(sizeof(curve__morton_to_hilbert3) == 96,
                    "morton_to_hilbert3 is 12 A4 states x 8 octants");
CURVE_STATIC_ASSERT(sizeof(curve__hilbert_to_morton3) == 96,
                    "hilbert_to_morton3 is 12 A4 states x 8 octants");

CURVE_INLINE uint32_t curve__transform3(uint32_t in, uint32_t bits,
                                        const uint8_t *CURVE_RESTRICT table)
{
    /* Table entries pack (next_state << 3) | output_digit: the low 3 bits are
     * the output octant, the upper bits the A4 state for the next level. */
    uint32_t transform = 0;
    uint32_t out = 0;
    int i;
    for (i = (int)(3 * (bits - 1)); i >= 0; i -= 3) {
        transform = table[transform | ((in >> i) & 7u)];
        out = (out << 3) | (transform & 7u);
        transform &= ~7u;
    }
    return out;
}

CURVE_INLINE uint32_t curve_hilbert3_xyz2d(uint32_t x, uint32_t y, uint32_t z, int bits)
{
    bits = curve__clamp_bits3(bits);
    x &= curve__lo_mask(bits);
    y &= curve__lo_mask(bits);
    z &= curve__lo_mask(bits);
    uint32_t m = curve__part1by2_u32(x) | (curve__part1by2_u32(y) << 1) |
                 (curve__part1by2_u32(z) << 2);
    return curve__transform3(m, (uint32_t)bits, curve__morton_to_hilbert3);
}

CURVE_INLINE void curve_hilbert3_d2xyz(uint32_t d, int bits,
                                       uint32_t *CURVE_RESTRICT x,
                                       uint32_t *CURVE_RESTRICT y,
                                       uint32_t *CURVE_RESTRICT z)
{
    bits = curve__clamp_bits3(bits);
    uint32_t m = curve__transform3(d, (uint32_t)bits, curve__hilbert_to_morton3);
    *x = curve__compact1by2_u32(m);
    *y = curve__compact1by2_u32(m >> 1);
    *z = curve__compact1by2_u32(m >> 2);
}

/* -------------------------------------------------------------------------- */
/* n-D Hilbert — Skilling's formulation of the Butz algorithm                 */
/*                                                                            */
/* coords[] holds n_dims coordinates of `bits` bits each; the Hilbert index   */
/* has n_dims*bits bits and is returned/passed as uint64_t, so require        */
/* n_dims*bits <= 64.  n_dims < 1 is a no-op (xy2d returns 0, d2xy writes      */
/* nothing).  n_dims is clamped down to CURVE_HILBERTND_MAX_DIMS and bits to   */
/* 1..min(32, 64/n_dims).  Coordinates are masked to `bits` bits.  Packing:   */
/* most significant level first; this matches curve_hilbert_xy2d(x, y) with   */
/* coords = {x, y}.  For n_dims=3 the orientation convention differs from     */
/* curve_hilbert3_xyz2d (both are valid Hilbert curves under different        */
/* rotations).                                                               */
/*                                                                            */
/* Reference: A. R. Butz, "Alternative Algorithm for Hilbert's Space-Filling  */
/* Curve", IEEE Trans. Comp. (1971); Doug Moore's Rice U. implementation;     */
/* J. Skilling, "Programming the Hilbert curve" (2004).                       */
/* -------------------------------------------------------------------------- */

#define CURVE_HILBERTND_MAX_DIMS 64

CURVE_INLINE void curve__hilbertnd_clamp(int *n_dims, int *bits)
{
    int maxbits;
    if (*n_dims < 1) *n_dims = 1;
    if (*n_dims > CURVE_HILBERTND_MAX_DIMS) *n_dims = CURVE_HILBERTND_MAX_DIMS;
    maxbits = 64 / *n_dims;
    if (maxbits > 32) maxbits = 32;
    if (*bits < 1) *bits = 1;
    if (*bits > maxbits) *bits = maxbits;
}

/* Skilling AxesToTranspose: in-place transform of X[0..n) (b bits each). */
CURVE_INLINE void curve__axes_to_transpose(uint32_t *X, int b, int n)
{
    uint32_t M = 1u << (b - 1), P, Q, t;
    int i;

    for (Q = M; Q > 1; Q >>= 1) {
        P = Q - 1;
        for (i = 0; i < n; i++) {
            if (X[i] & Q) {
                X[0] ^= P;
            } else {
                t = (X[0] ^ X[i]) & P;
                X[0] ^= t;
                X[i] ^= t;
            }
        }
    }
    /* Gray encode */
    for (i = 1; i < n; i++) X[i] ^= X[i - 1];
    t = 0;
    for (Q = M; Q > 1; Q >>= 1)
        if (X[n - 1] & Q) t ^= Q - 1;
    for (i = 0; i < n; i++) X[i] ^= t;
}

/* Skilling TransposeToAxes: inverse of curve__axes_to_transpose. */
CURVE_INLINE void curve__transpose_to_axes(uint32_t *X, int b, int n)
{
    uint32_t P, Q, t;
    int i, k;

    /* Gray decode by H ^ (H/2) */
    t = X[n - 1] >> 1;
    for (i = n - 1; i > 0; i--) X[i] ^= X[i - 1];
    X[0] ^= t;

    /* Q = 2, 4, ..., 2^(b-1).  Counted so b == 32 does not need 2^b in a
     * uint32_t (the old `Q != (2u << (b-1))` relied on unsigned wrap to 0). */
    Q = 2u;
    for (k = 1; k < b; k++) {
        P = Q - 1u;
        for (i = n - 1; i >= 0; i--) {
            if (X[i] & Q) {
                X[0] ^= P;
            } else {
                t = (X[0] ^ X[i]) & P;
                X[0] ^= t;
                X[i] ^= t;
            }
        }
        Q <<= 1;
    }
}

CURVE_INLINE uint64_t curve_hilbertnd_xy2d(const uint32_t *CURVE_RESTRICT coords,
                                           int n_dims, int bits)
{
    uint32_t X[CURVE_HILBERTND_MAX_DIMS];
    uint64_t h = 0;
    int i, p;
    uint32_t mask;

    /* Do not bump n_dims upward: n_dims < 1 would otherwise read coords[0]. */
    if (n_dims < 1) return 0;
    curve__hilbertnd_clamp(&n_dims, &bits);
    mask = curve__lo_mask(bits);
    for (i = 0; i < n_dims; i++) X[i] = coords[i] & mask;
    if (n_dims > 1)
        curve__axes_to_transpose(X, bits, n_dims);

    for (p = bits - 1; p >= 0; p--)
        for (i = 0; i < n_dims; i++)
            h = (h << 1) | ((X[i] >> p) & 1u);
    return h;
}

CURVE_INLINE void curve_hilbertnd_d2xy(uint64_t d, uint32_t *CURVE_RESTRICT coords,
                                       int n_dims, int bits)
{
    uint32_t X[CURVE_HILBERTND_MAX_DIMS];
    int i, p;

    if (n_dims < 1) return;
    curve__hilbertnd_clamp(&n_dims, &bits);
    for (i = 0; i < n_dims; i++) X[i] = 0;

    for (p = bits - 1; p >= 0; p--) {
        for (i = 0; i < n_dims; i++) {
            X[i] = (X[i] << 1) | (uint32_t)((d >> (n_dims * bits - 1)) & 1u);
            d <<= 1;
        }
    }
    /* X[i] now holds the transpose word of dim i (dim 0 = MSB of each level) */
    if (n_dims > 1)
        curve__transpose_to_axes(X, bits, n_dims);
    for (i = 0; i < n_dims; i++) coords[i] = X[i];
}

/* -------------------------------------------------------------------------- */
/* Sierpiński curve ordering key (Bartholdi & Platzman 1989), bits in 1..31   */
/*                                                                            */
/* Maps grid point (x, y), with x,y < 2^bits, to an ordering key in           */
/* [0, 2^(2*bits+1)).  Sorting points by key yields Sierpiński curve order.   */
/* The key is NOT dense (half the key space is unused: the curve is closed    */
/* and its two triangular halves share the main diagonal), so no inverse      */
/* d2xy is provided.  Keys are injective over the grid.                       */
/* -------------------------------------------------------------------------- */

CURVE_INLINE int curve__clamp_bitss(int bits)
{
    if (bits < 1)  return 1;
    if (bits > 31) return 31;
    return bits;
}

CURVE_INLINE uint64_t curve_sierpinski_xy2d(uint32_t x, uint32_t y, int bits)
{
    uint64_t result = 0;
    uint32_t m, loop, t;

    bits = curve__clamp_bitss(bits);
    x &= curve__lo_mask(bits);
    y &= curve__lo_mask(bits);
    m = (1u << bits) - 1u;

    if (x > y) {
        result = 1;
        x = m - x;   /* point-reflect through center (Bartholdi) */
        y = m - y;
    }
    for (loop = m; loop > 0; loop >>= 1) {
        result <<= 1;
        if (x + y > m) {
            result |= 1;
            t = x; x = m - y; y = t;
        }
        x <<= 1;
        y <<= 1;
        result <<= 1;
        if (y > m) {
            result |= 1;
            t = x; x = y - m; y = m - t;
        }
    }
    return result;
}

/* -------------------------------------------------------------------------- */
/* Gilbert curves — Hamiltonian paths on arbitrary rectangles                 */
/*                                                                            */
/* Exact-cover visitors for arbitrary (non-power-of-two) rectangles/boxes.    */
/* 2D: consecutive cells are always edge-adjacent.  3D: consecutive cells are */
/* face-adjacent when all box dimensions are even; with odd dimensions the    */
/* joins between sub-blocks may contain diagonal steps (as in the reference   */
/* implementation — "even sizes are recommended in 3D").                      */
/* -------------------------------------------------------------------------- */

typedef void (*curve_visit2_fn)(int x, int y, void *user);
typedef void (*curve_visit3_fn)(int x, int y, int z, void *user);

CURVE_INLINE int curve__sgn(int x)
{
    return (x < 0) ? -1 : (x > 0) ? 1 : 0;
}

CURVE_INLINE void curve__gilbert2_r(int x, int y, int ax, int ay, int bx, int by,
                                    curve_visit2_fn cb, void *user)
{
    int w = ax + ay; if (w < 0) w = -w;
    int h = bx + by; if (h < 0) h = -h;
    int dax = curve__sgn(ax), day = curve__sgn(ay);
    int dbx = curve__sgn(bx), dby = curve__sgn(by);
    int i;

    if (h == 1) {
        for (i = 0; i < w; ++i) {
            cb(x, y, user);
            x += dax; y += day;
        }
        return;
    }
    if (w == 1) {
        for (i = 0; i < h; ++i) {
            cb(x, y, user);
            x += dbx; y += dby;
        }
        return;
    }

    int ax2 = ax / 2, ay2 = ay / 2;
    int bx2 = bx / 2, by2 = by / 2;
    int w2 = ax2 + ay2; if (w2 < 0) w2 = -w2;
    int h2 = bx2 + by2; if (h2 < 0) h2 = -h2;

    if (2 * w > 3 * h) {
        if ((w2 & 1) && w > 2) { ax2 += dax; ay2 += day; }
        curve__gilbert2_r(x, y, ax2, ay2, bx, by, cb, user);
        curve__gilbert2_r(x + ax2, y + ay2, ax - ax2, ay - ay2, bx, by, cb, user);
    } else {
        if ((h2 & 1) && h > 2) { bx2 += dbx; by2 += dby; }
        curve__gilbert2_r(x, y, bx2, by2, ax2, ay2, cb, user);
        curve__gilbert2_r(x + bx2, y + by2, ax, ay, bx - bx2, by - by2, cb, user);
        curve__gilbert2_r(x + (ax - dax) + (bx2 - dbx),
                          y + (ay - day) + (by2 - dby),
                          -bx2, -by2, -(ax - ax2), -(ay - ay2), cb, user);
    }
}

CURVE_INLINE void curve_gilbert2(int width, int height, curve_visit2_fn cb, void *user)
{
    if (CURVE_UNLIKELY(!cb || width <= 0 || height <= 0)) return;
    if (width >= height)
        curve__gilbert2_r(0, 0, width, 0, 0, height, cb, user);
    else
        curve__gilbert2_r(0, 0, 0, height, width, 0, cb, user);
}

CURVE_INLINE void curve__gilbert3_r(int x, int y, int z,
                                    int ax, int ay, int az,
                                    int bx, int by, int bz,
                                    int cx, int cy, int cz,
                                    curve_visit3_fn cb, void *user)
{
    int w = ax + ay + az; if (w < 0) w = -w;
    int h = bx + by + bz; if (h < 0) h = -h;
    int d = cx + cy + cz; if (d < 0) d = -d;

    int dax = curve__sgn(ax), day = curve__sgn(ay), daz = curve__sgn(az);
    int dbx = curve__sgn(bx), dby = curve__sgn(by), dbz = curve__sgn(bz);
    int dcx = curve__sgn(cx), dcy = curve__sgn(cy), dcz = curve__sgn(cz);
    int i;

    if (h == 1 && d == 1) {
        for (i = 0; i < w; ++i) { cb(x, y, z, user); x += dax; y += day; z += daz; }
        return;
    }
    if (w == 1 && d == 1) {
        for (i = 0; i < h; ++i) { cb(x, y, z, user); x += dbx; y += dby; z += dbz; }
        return;
    }
    if (w == 1 && h == 1) {
        for (i = 0; i < d; ++i) { cb(x, y, z, user); x += dcx; y += dcy; z += dcz; }
        return;
    }

    int ax2 = ax / 2, ay2 = ay / 2, az2 = az / 2;
    int bx2 = bx / 2, by2 = by / 2, bz2 = bz / 2;
    int cx2 = cx / 2, cy2 = cy / 2, cz2 = cz / 2;

    int w2 = ax2 + ay2 + az2; if (w2 < 0) w2 = -w2;
    int h2 = bx2 + by2 + bz2; if (h2 < 0) h2 = -h2;
    int d2 = cx2 + cy2 + cz2; if (d2 < 0) d2 = -d2;

    if ((w2 & 1) && w > 2) { ax2 += dax; ay2 += day; az2 += daz; }
    if ((h2 & 1) && h > 2) { bx2 += dbx; by2 += dby; bz2 += dbz; }
    if ((d2 & 1) && d > 2) { cx2 += dcx; cy2 += dcy; cz2 += dcz; }

    if ((2 * w > 3 * h) && (2 * w > 3 * d)) {
        curve__gilbert3_r(x, y, z, ax2, ay2, az2, bx, by, bz, cx, cy, cz, cb, user);
        curve__gilbert3_r(x + ax2, y + ay2, z + az2,
                          ax - ax2, ay - ay2, az - az2, bx, by, bz, cx, cy, cz, cb, user);
    } else if (3 * h > 4 * d) {
        curve__gilbert3_r(x, y, z, bx2, by2, bz2, cx, cy, cz, ax2, ay2, az2, cb, user);
        curve__gilbert3_r(x + bx2, y + by2, z + bz2,
                          ax, ay, az, bx - bx2, by - by2, bz - bz2, cx, cy, cz, cb, user);
        curve__gilbert3_r(x + (ax - dax) + (bx2 - dbx),
                          y + (ay - day) + (by2 - dby),
                          z + (az - daz) + (bz2 - dbz),
                          -bx2, -by2, -bz2, cx, cy, cz,
                          -(ax - ax2), -(ay - ay2), -(az - az2), cb, user);
    } else if (3 * d > 4 * h) {
        curve__gilbert3_r(x, y, z, cx2, cy2, cz2, ax2, ay2, az2, bx, by, bz, cb, user);
        curve__gilbert3_r(x + cx2, y + cy2, z + cz2,
                          ax, ay, az, bx, by, bz, cx - cx2, cy - cy2, cz - cz2, cb, user);
        curve__gilbert3_r(x + (ax - dax) + (cx2 - dcx),
                          y + (ay - day) + (cy2 - dcy),
                          z + (az - daz) + (cz2 - dcz),
                          -cx2, -cy2, -cz2, -(ax - ax2), -(ay - ay2), -(az - az2),
                          bx, by, bz, cb, user);
    } else {
        curve__gilbert3_r(x, y, z, bx2, by2, bz2, cx2, cy2, cz2, ax2, ay2, az2, cb, user);
        curve__gilbert3_r(x + bx2, y + by2, z + bz2,
                          cx, cy, cz, ax2, ay2, az2, bx - bx2, by - by2, bz - bz2, cb, user);
        curve__gilbert3_r(x + (bx2 - dbx) + (cx - dcx),
                          y + (by2 - dby) + (cy - dcy),
                          z + (bz2 - dbz) + (cz - dcz),
                          ax, ay, az, -bx2, -by2, -bz2,
                          -(cx - cx2), -(cy - cy2), -(cz - cz2), cb, user);
        curve__gilbert3_r(x + (ax - dax) + bx2 + (cx - dcx),
                          y + (ay - day) + by2 + (cy - dcy),
                          z + (az - daz) + bz2 + (cz - dcz),
                          -cx, -cy, -cz, -(ax - ax2), -(ay - ay2), -(az - az2),
                          bx - bx2, by - by2, bz - bz2, cb, user);
        curve__gilbert3_r(x + (ax - dax) + (bx2 - dbx),
                          y + (ay - day) + (by2 - dby),
                          z + (az - daz) + (bz2 - dbz),
                          -bx2, -by2, -bz2, cx2, cy2, cz2,
                          -(ax - ax2), -(ay - ay2), -(az - az2), cb, user);
    }
}

CURVE_INLINE void curve_gilbert3(int width, int height, int depth,
                                 curve_visit3_fn cb, void *user)
{
    if (CURVE_UNLIKELY(!cb || width <= 0 || height <= 0 || depth <= 0)) return;
    if (width >= height && width >= depth)
        curve__gilbert3_r(0, 0, 0, width, 0, 0, 0, height, 0, 0, 0, depth, cb, user);
    else if (height >= width && height >= depth)
        curve__gilbert3_r(0, 0, 0, 0, height, 0, width, 0, 0, 0, 0, depth, cb, user);
    else
        curve__gilbert3_r(0, 0, 0, 0, 0, depth, width, 0, 0, 0, height, 0, cb, user);
}

/* -------------------------------------------------------------------------- */
/* SIMD / batch APIs                                                          */
/*                                                                            */
/* SoA arrays, no aliasing (restrict).  Unaligned loads/stores.  Widest ISA   */
/* present is used; a scalar tail handles leftovers.  32-bit Morton 2D uses   */
/* NEON / AVX2 / AVX-512; 64-bit Morton 2D/3D batch uses AVX2 / AVX-512 only  */
/* (2-wide NEON SWAR loses to the compiler's scalar loop on aarch64).  Scalar */
/* BMI2 PDEP/PEXT still applies to the single-element APIs above.             */
/* -------------------------------------------------------------------------- */

CURVE_INLINE const char *curve_simd_isa(void)
{
#if defined(CURVE_HAVE_AVX512) && defined(CURVE_HAVE_BMI2)
    return "avx512+bmi2";
#elif defined(CURVE_HAVE_AVX512)
    return "avx512";
#elif defined(CURVE_HAVE_AVX2) && defined(CURVE_HAVE_BMI2)
    return "avx2+bmi2";
#elif defined(CURVE_HAVE_AVX2)
    return "avx2";
#elif defined(CURVE_HAVE_BMI2)
    return "bmi2";
#elif defined(CURVE_HAVE_NEON)
    return "neon";
#else
    return "scalar";
#endif
}

#if defined(CURVE_HAVE_NEON)
CURVE_INLINE uint32x4_t curve__part1by1_u32x4(uint32x4_t x)
{
    x = vandq_u32(x, vdupq_n_u32(0x0000FFFFu));
    x = vandq_u32(vorrq_u32(x, vshlq_n_u32(x, 8)), vdupq_n_u32(0x00FF00FFu));
    x = vandq_u32(vorrq_u32(x, vshlq_n_u32(x, 4)), vdupq_n_u32(0x0F0F0F0Fu));
    x = vandq_u32(vorrq_u32(x, vshlq_n_u32(x, 2)), vdupq_n_u32(0x33333333u));
    x = vandq_u32(vorrq_u32(x, vshlq_n_u32(x, 1)), vdupq_n_u32(0x55555555u));
    return x;
}
CURVE_INLINE uint32x4_t curve__compact1by1_u32x4(uint32x4_t x)
{
    x = vandq_u32(x, vdupq_n_u32(0x55555555u));
    x = vandq_u32(vorrq_u32(x, vshrq_n_u32(x, 1)), vdupq_n_u32(0x33333333u));
    x = vandq_u32(vorrq_u32(x, vshrq_n_u32(x, 2)), vdupq_n_u32(0x0F0F0F0Fu));
    x = vandq_u32(vorrq_u32(x, vshrq_n_u32(x, 4)), vdupq_n_u32(0x00FF00FFu));
    x = vandq_u32(vorrq_u32(x, vshrq_n_u32(x, 8)), vdupq_n_u32(0x0000FFFFu));
    return x;
}
#endif /* CURVE_HAVE_NEON */

#if defined(CURVE_HAVE_AVX2)
CURVE_INLINE __m256i curve__part1by1_avx2(__m256i x)
{
    x = _mm256_and_si256(x, _mm256_set1_epi32((int)0x0000FFFFu));
    x = _mm256_and_si256(_mm256_or_si256(x, _mm256_slli_epi32(x, 8)),
                         _mm256_set1_epi32((int)0x00FF00FFu));
    x = _mm256_and_si256(_mm256_or_si256(x, _mm256_slli_epi32(x, 4)),
                         _mm256_set1_epi32((int)0x0F0F0F0Fu));
    x = _mm256_and_si256(_mm256_or_si256(x, _mm256_slli_epi32(x, 2)),
                         _mm256_set1_epi32((int)0x33333333u));
    x = _mm256_and_si256(_mm256_or_si256(x, _mm256_slli_epi32(x, 1)),
                         _mm256_set1_epi32((int)0x55555555u));
    return x;
}
CURVE_INLINE __m256i curve__compact1by1_avx2(__m256i x)
{
    x = _mm256_and_si256(x, _mm256_set1_epi32((int)0x55555555u));
    x = _mm256_and_si256(_mm256_or_si256(x, _mm256_srli_epi32(x, 1)),
                         _mm256_set1_epi32((int)0x33333333u));
    x = _mm256_and_si256(_mm256_or_si256(x, _mm256_srli_epi32(x, 2)),
                         _mm256_set1_epi32((int)0x0F0F0F0Fu));
    x = _mm256_and_si256(_mm256_or_si256(x, _mm256_srli_epi32(x, 4)),
                         _mm256_set1_epi32((int)0x00FF00FFu));
    x = _mm256_and_si256(_mm256_or_si256(x, _mm256_srli_epi32(x, 8)),
                         _mm256_set1_epi32((int)0x0000FFFFu));
    return x;
}
CURVE_INLINE __m256i curve__part1by1_avx2_64(__m256i x)
{
    x = _mm256_and_si256(x, _mm256_set1_epi64x((long long)0x00000000ffffffffull));
    x = _mm256_and_si256(_mm256_or_si256(x, _mm256_slli_epi64(x, 16)),
                         _mm256_set1_epi64x((long long)0x0000ffff0000ffffull));
    x = _mm256_and_si256(_mm256_or_si256(x, _mm256_slli_epi64(x, 8)),
                         _mm256_set1_epi64x((long long)0x00ff00ff00ff00ffull));
    x = _mm256_and_si256(_mm256_or_si256(x, _mm256_slli_epi64(x, 4)),
                         _mm256_set1_epi64x((long long)0x0f0f0f0f0f0f0f0full));
    x = _mm256_and_si256(_mm256_or_si256(x, _mm256_slli_epi64(x, 2)),
                         _mm256_set1_epi64x((long long)0x3333333333333333ull));
    x = _mm256_and_si256(_mm256_or_si256(x, _mm256_slli_epi64(x, 1)),
                         _mm256_set1_epi64x((long long)0x5555555555555555ull));
    return x;
}
CURVE_INLINE __m256i curve__compact1by1_avx2_64(__m256i x)
{
    x = _mm256_and_si256(x, _mm256_set1_epi64x((long long)0x5555555555555555ull));
    x = _mm256_and_si256(_mm256_or_si256(x, _mm256_srli_epi64(x, 1)),
                         _mm256_set1_epi64x((long long)0x3333333333333333ull));
    x = _mm256_and_si256(_mm256_or_si256(x, _mm256_srli_epi64(x, 2)),
                         _mm256_set1_epi64x((long long)0x0f0f0f0f0f0f0f0full));
    x = _mm256_and_si256(_mm256_or_si256(x, _mm256_srli_epi64(x, 4)),
                         _mm256_set1_epi64x((long long)0x00ff00ff00ff00ffull));
    x = _mm256_and_si256(_mm256_or_si256(x, _mm256_srli_epi64(x, 8)),
                         _mm256_set1_epi64x((long long)0x0000ffff0000ffffull));
    x = _mm256_and_si256(_mm256_or_si256(x, _mm256_srli_epi64(x, 16)),
                         _mm256_set1_epi64x((long long)0x00000000ffffffffull));
    return x;
}
CURVE_INLINE __m256i curve__part1by2_avx2_64(__m256i x)
{
    x = _mm256_and_si256(x, _mm256_set1_epi64x((long long)0x1fffffull));
    x = _mm256_and_si256(_mm256_or_si256(x, _mm256_slli_epi64(x, 32)),
                         _mm256_set1_epi64x((long long)0x1f00000000ffffull));
    x = _mm256_and_si256(_mm256_or_si256(x, _mm256_slli_epi64(x, 16)),
                         _mm256_set1_epi64x((long long)0x1f0000ff0000ffull));
    x = _mm256_and_si256(_mm256_or_si256(x, _mm256_slli_epi64(x, 8)),
                         _mm256_set1_epi64x((long long)0x100f00f00f00f00full));
    x = _mm256_and_si256(_mm256_or_si256(x, _mm256_slli_epi64(x, 4)),
                         _mm256_set1_epi64x((long long)0x10c30c30c30c30c3ull));
    x = _mm256_and_si256(_mm256_or_si256(x, _mm256_slli_epi64(x, 2)),
                         _mm256_set1_epi64x((long long)0x1249249249249249ull));
    return x;
}
CURVE_INLINE __m256i curve__compact1by2_avx2_64(__m256i x)
{
    x = _mm256_and_si256(x, _mm256_set1_epi64x((long long)0x1249249249249249ull));
    x = _mm256_and_si256(_mm256_or_si256(x, _mm256_srli_epi64(x, 2)),
                         _mm256_set1_epi64x((long long)0x10c30c30c30c30c3ull));
    x = _mm256_and_si256(_mm256_or_si256(x, _mm256_srli_epi64(x, 4)),
                         _mm256_set1_epi64x((long long)0x100f00f00f00f00full));
    x = _mm256_and_si256(_mm256_or_si256(x, _mm256_srli_epi64(x, 8)),
                         _mm256_set1_epi64x((long long)0x1f0000ff0000ffull));
    x = _mm256_and_si256(_mm256_or_si256(x, _mm256_srli_epi64(x, 16)),
                         _mm256_set1_epi64x((long long)0x1f00000000ffffull));
    x = _mm256_and_si256(_mm256_or_si256(x, _mm256_srli_epi64(x, 32)),
                         _mm256_set1_epi64x((long long)0x1fffffull));
    return x;
}
CURVE_INLINE __m128i curve__pack64lo32_avx2(__m256i v)
{
    /* Low 32 bits of each 64-bit lane -> 4x uint32. */
    v = _mm256_shuffle_epi32(v, _MM_SHUFFLE(2, 0, 2, 0));
    return _mm_unpacklo_epi64(_mm256_castsi256_si128(v),
                              _mm256_extracti128_si256(v, 1));
}
#endif /* CURVE_HAVE_AVX2 */

#if defined(CURVE_HAVE_AVX512)
CURVE_INLINE __m512i curve__part1by1_avx512(__m512i x)
{
    x = _mm512_and_si512(x, _mm512_set1_epi32((int)0x0000FFFFu));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_slli_epi32(x, 8)),
                         _mm512_set1_epi32((int)0x00FF00FFu));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_slli_epi32(x, 4)),
                         _mm512_set1_epi32((int)0x0F0F0F0Fu));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_slli_epi32(x, 2)),
                         _mm512_set1_epi32((int)0x33333333u));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_slli_epi32(x, 1)),
                         _mm512_set1_epi32((int)0x55555555u));
    return x;
}
CURVE_INLINE __m512i curve__compact1by1_avx512(__m512i x)
{
    x = _mm512_and_si512(x, _mm512_set1_epi32((int)0x55555555u));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_srli_epi32(x, 1)),
                         _mm512_set1_epi32((int)0x33333333u));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_srli_epi32(x, 2)),
                         _mm512_set1_epi32((int)0x0F0F0F0Fu));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_srli_epi32(x, 4)),
                         _mm512_set1_epi32((int)0x00FF00FFu));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_srli_epi32(x, 8)),
                         _mm512_set1_epi32((int)0x0000FFFFu));
    return x;
}
CURVE_INLINE __m512i curve__part1by1_avx512_64(__m512i x)
{
    x = _mm512_and_si512(x, _mm512_set1_epi64((long long)0x00000000ffffffffull));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_slli_epi64(x, 16)),
                         _mm512_set1_epi64((long long)0x0000ffff0000ffffull));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_slli_epi64(x, 8)),
                         _mm512_set1_epi64((long long)0x00ff00ff00ff00ffull));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_slli_epi64(x, 4)),
                         _mm512_set1_epi64((long long)0x0f0f0f0f0f0f0f0full));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_slli_epi64(x, 2)),
                         _mm512_set1_epi64((long long)0x3333333333333333ull));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_slli_epi64(x, 1)),
                         _mm512_set1_epi64((long long)0x5555555555555555ull));
    return x;
}
CURVE_INLINE __m512i curve__compact1by1_avx512_64(__m512i x)
{
    x = _mm512_and_si512(x, _mm512_set1_epi64((long long)0x5555555555555555ull));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_srli_epi64(x, 1)),
                         _mm512_set1_epi64((long long)0x3333333333333333ull));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_srli_epi64(x, 2)),
                         _mm512_set1_epi64((long long)0x0f0f0f0f0f0f0f0full));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_srli_epi64(x, 4)),
                         _mm512_set1_epi64((long long)0x00ff00ff00ff00ffull));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_srli_epi64(x, 8)),
                         _mm512_set1_epi64((long long)0x0000ffff0000ffffull));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_srli_epi64(x, 16)),
                         _mm512_set1_epi64((long long)0x00000000ffffffffull));
    return x;
}
CURVE_INLINE __m512i curve__part1by2_avx512_64(__m512i x)
{
    x = _mm512_and_si512(x, _mm512_set1_epi64((long long)0x1fffffull));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_slli_epi64(x, 32)),
                         _mm512_set1_epi64((long long)0x1f00000000ffffull));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_slli_epi64(x, 16)),
                         _mm512_set1_epi64((long long)0x1f0000ff0000ffull));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_slli_epi64(x, 8)),
                         _mm512_set1_epi64((long long)0x100f00f00f00f00full));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_slli_epi64(x, 4)),
                         _mm512_set1_epi64((long long)0x10c30c30c30c30c3ull));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_slli_epi64(x, 2)),
                         _mm512_set1_epi64((long long)0x1249249249249249ull));
    return x;
}
CURVE_INLINE __m512i curve__compact1by2_avx512_64(__m512i x)
{
    x = _mm512_and_si512(x, _mm512_set1_epi64((long long)0x1249249249249249ull));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_srli_epi64(x, 2)),
                         _mm512_set1_epi64((long long)0x10c30c30c30c30c3ull));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_srli_epi64(x, 4)),
                         _mm512_set1_epi64((long long)0x100f00f00f00f00full));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_srli_epi64(x, 8)),
                         _mm512_set1_epi64((long long)0x1f0000ff0000ffull));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_srli_epi64(x, 16)),
                         _mm512_set1_epi64((long long)0x1f00000000ffffull));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_srli_epi64(x, 32)),
                         _mm512_set1_epi64((long long)0x1fffffull));
    return x;
}
#endif /* CURVE_HAVE_AVX512 */

CURVE_INLINE void curve_gray_encode_batch(const uint64_t *CURVE_RESTRICT x,
                                          uint64_t *CURVE_RESTRICT out, size_t n)
{
    size_t i = 0;
#if defined(CURVE_HAVE_AVX512)
    for (; i < n - (n & 7u); i += 8) {
        __m512i v = _mm512_loadu_si512((const void *)(x + i));
        _mm512_storeu_si512((void *)(out + i),
                            _mm512_xor_si512(v, _mm512_srli_epi64(v, 1)));
    }
#elif defined(CURVE_HAVE_AVX2)
    for (; i < n - (n & 3u); i += 4) {
        __m256i v = _mm256_loadu_si256((const __m256i *)(const void *)(x + i));
        _mm256_storeu_si256((__m256i *)(void *)(out + i),
                            _mm256_xor_si256(v, _mm256_srli_epi64(v, 1)));
    }
#elif defined(CURVE_HAVE_NEON)
    for (; i < n - (n & 1u); i += 2) {
        uint64x2_t v = vld1q_u64(x + i);
        vst1q_u64(out + i, veorq_u64(v, vshrq_n_u64(v, 1)));
    }
#endif
    for (; i < n; ++i) out[i] = curve_gray_encode(x[i]);
}

CURVE_INLINE void curve_gray_decode_batch(const uint64_t *CURVE_RESTRICT g,
                                          uint64_t *CURVE_RESTRICT out, size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) out[i] = curve_gray_decode(g[i]);
}

CURVE_INLINE void curve_morton2_batch(const uint32_t *CURVE_RESTRICT x,
                                      const uint32_t *CURVE_RESTRICT y,
                                      uint32_t *CURVE_RESTRICT out, size_t n)
{
    size_t i = 0;
#if defined(CURVE_HAVE_AVX512)
    for (; i < n - (n & 15u); i += 16) {
        __m512i vx = _mm512_loadu_si512((const void *)(x + i));
        __m512i vy = _mm512_loadu_si512((const void *)(y + i));
        _mm512_storeu_si512((void *)(out + i),
            _mm512_or_si512(curve__part1by1_avx512(vx),
                            _mm512_slli_epi32(curve__part1by1_avx512(vy), 1)));
    }
#elif defined(CURVE_HAVE_AVX2)
    for (; i < n - (n & 7u); i += 8) {
        __m256i vx = _mm256_loadu_si256((const __m256i *)(const void *)(x + i));
        __m256i vy = _mm256_loadu_si256((const __m256i *)(const void *)(y + i));
        _mm256_storeu_si256((__m256i *)(void *)(out + i),
            _mm256_or_si256(curve__part1by1_avx2(vx),
                            _mm256_slli_epi32(curve__part1by1_avx2(vy), 1)));
    }
#elif defined(CURVE_HAVE_NEON)
    for (; i < n - (n & 3u); i += 4) {
        uint32x4_t vx = vld1q_u32(x + i);
        uint32x4_t vy = vld1q_u32(y + i);
        vst1q_u32(out + i, vorrq_u32(curve__part1by1_u32x4(vx),
                                     vshlq_n_u32(curve__part1by1_u32x4(vy), 1)));
    }
#endif
    for (; i < n; ++i) out[i] = curve_morton2(x[i], y[i]);
}

CURVE_INLINE void curve_morton2_decode_batch(const uint32_t *CURVE_RESTRICT d,
                                             uint32_t *CURVE_RESTRICT x,
                                             uint32_t *CURVE_RESTRICT y,
                                             size_t n)
{
    size_t i = 0;
#if defined(CURVE_HAVE_AVX512)
    for (; i < n - (n & 15u); i += 16) {
        __m512i vd = _mm512_loadu_si512((const void *)(d + i));
        _mm512_storeu_si512((void *)(x + i), curve__compact1by1_avx512(vd));
        _mm512_storeu_si512((void *)(y + i),
                            curve__compact1by1_avx512(_mm512_srli_epi32(vd, 1)));
    }
#elif defined(CURVE_HAVE_AVX2)
    for (; i < n - (n & 7u); i += 8) {
        __m256i vd = _mm256_loadu_si256((const __m256i *)(const void *)(d + i));
        _mm256_storeu_si256((__m256i *)(void *)(x + i), curve__compact1by1_avx2(vd));
        _mm256_storeu_si256((__m256i *)(void *)(y + i),
                            curve__compact1by1_avx2(_mm256_srli_epi32(vd, 1)));
    }
#elif defined(CURVE_HAVE_NEON)
    for (; i < n - (n & 3u); i += 4) {
        uint32x4_t vd = vld1q_u32(d + i);
        vst1q_u32(x + i, curve__compact1by1_u32x4(vd));
        vst1q_u32(y + i, curve__compact1by1_u32x4(vshrq_n_u32(vd, 1)));
    }
#endif
    for (; i < n; ++i) curve_morton2_decode(d[i], x + i, y + i);
}

CURVE_INLINE void curve_morton2_64_batch(const uint32_t *CURVE_RESTRICT x,
                                         const uint32_t *CURVE_RESTRICT y,
                                         uint64_t *CURVE_RESTRICT out, size_t n)
{
    size_t i = 0;
#if defined(CURVE_HAVE_AVX512)
    for (; i < n - (n & 7u); i += 8) {
        __m256i x32 = _mm256_loadu_si256((const __m256i *)(const void *)(x + i));
        __m256i y32 = _mm256_loadu_si256((const __m256i *)(const void *)(y + i));
        __m512i vx = curve__part1by1_avx512_64(_mm512_cvtepu32_epi64(x32));
        __m512i vy = curve__part1by1_avx512_64(_mm512_cvtepu32_epi64(y32));
        _mm512_storeu_si512((void *)(out + i),
                            _mm512_or_si512(vx, _mm512_slli_epi64(vy, 1)));
    }
#elif defined(CURVE_HAVE_AVX2)
    for (; i < n - (n & 3u); i += 4) {
        __m128i x32 = _mm_loadu_si128((const __m128i *)(const void *)(x + i));
        __m128i y32 = _mm_loadu_si128((const __m128i *)(const void *)(y + i));
        __m256i vx = curve__part1by1_avx2_64(_mm256_cvtepu32_epi64(x32));
        __m256i vy = curve__part1by1_avx2_64(_mm256_cvtepu32_epi64(y32));
        _mm256_storeu_si256((__m256i *)(void *)(out + i),
                            _mm256_or_si256(vx, _mm256_slli_epi64(vy, 1)));
    }
#endif
    for (; i < n; ++i) out[i] = curve_morton2_64(x[i], y[i]);
}

CURVE_INLINE void curve_morton2_64_decode_batch(const uint64_t *CURVE_RESTRICT d,
                                                uint32_t *CURVE_RESTRICT x,
                                                uint32_t *CURVE_RESTRICT y,
                                                size_t n)
{
    size_t i = 0;
#if defined(CURVE_HAVE_AVX512)
    for (; i < n - (n & 7u); i += 8) {
        __m512i vd = _mm512_loadu_si512((const void *)(d + i));
        _mm256_storeu_si256((__m256i *)(void *)(x + i),
            _mm512_cvtepi64_epi32(curve__compact1by1_avx512_64(vd)));
        _mm256_storeu_si256((__m256i *)(void *)(y + i),
            _mm512_cvtepi64_epi32(curve__compact1by1_avx512_64(_mm512_srli_epi64(vd, 1))));
    }
#elif defined(CURVE_HAVE_AVX2)
    for (; i < n - (n & 3u); i += 4) {
        __m256i vd = _mm256_loadu_si256((const __m256i *)(const void *)(d + i));
        _mm_storeu_si128((__m128i *)(void *)(x + i),
            curve__pack64lo32_avx2(curve__compact1by1_avx2_64(vd)));
        _mm_storeu_si128((__m128i *)(void *)(y + i),
            curve__pack64lo32_avx2(curve__compact1by1_avx2_64(_mm256_srli_epi64(vd, 1))));
    }
#endif
    for (; i < n; ++i) curve_morton2_64_decode(d[i], x + i, y + i);
}

CURVE_INLINE void curve_morton3_batch(const uint32_t *CURVE_RESTRICT x,
                                      const uint32_t *CURVE_RESTRICT y,
                                      const uint32_t *CURVE_RESTRICT z,
                                      uint64_t *CURVE_RESTRICT out, size_t n)
{
    size_t i = 0;
#if defined(CURVE_HAVE_AVX512)
    for (; i < n - (n & 7u); i += 8) {
        __m512i vx = curve__part1by2_avx512_64(
            _mm512_cvtepu32_epi64(_mm256_loadu_si256((const __m256i *)(const void *)(x + i))));
        __m512i vy = curve__part1by2_avx512_64(
            _mm512_cvtepu32_epi64(_mm256_loadu_si256((const __m256i *)(const void *)(y + i))));
        __m512i vz = curve__part1by2_avx512_64(
            _mm512_cvtepu32_epi64(_mm256_loadu_si256((const __m256i *)(const void *)(z + i))));
        _mm512_storeu_si512((void *)(out + i),
            _mm512_or_si512(vx, _mm512_or_si512(_mm512_slli_epi64(vy, 1),
                                                _mm512_slli_epi64(vz, 2))));
    }
#elif defined(CURVE_HAVE_AVX2)
    for (; i < n - (n & 3u); i += 4) {
        __m256i vx = curve__part1by2_avx2_64(
            _mm256_cvtepu32_epi64(_mm_loadu_si128((const __m128i *)(const void *)(x + i))));
        __m256i vy = curve__part1by2_avx2_64(
            _mm256_cvtepu32_epi64(_mm_loadu_si128((const __m128i *)(const void *)(y + i))));
        __m256i vz = curve__part1by2_avx2_64(
            _mm256_cvtepu32_epi64(_mm_loadu_si128((const __m128i *)(const void *)(z + i))));
        _mm256_storeu_si256((__m256i *)(void *)(out + i),
            _mm256_or_si256(vx, _mm256_or_si256(_mm256_slli_epi64(vy, 1),
                                                _mm256_slli_epi64(vz, 2))));
    }
#endif
    for (; i < n; ++i) out[i] = curve_morton3(x[i], y[i], z[i]);
}

CURVE_INLINE void curve_morton3_decode_batch(const uint64_t *CURVE_RESTRICT d,
                                             uint32_t *CURVE_RESTRICT x,
                                             uint32_t *CURVE_RESTRICT y,
                                             uint32_t *CURVE_RESTRICT z,
                                             size_t n)
{
    size_t i = 0;
#if defined(CURVE_HAVE_AVX512)
    for (; i < n - (n & 7u); i += 8) {
        __m512i vd = _mm512_loadu_si512((const void *)(d + i));
        _mm256_storeu_si256((__m256i *)(void *)(x + i),
            _mm512_cvtepi64_epi32(curve__compact1by2_avx512_64(vd)));
        _mm256_storeu_si256((__m256i *)(void *)(y + i),
            _mm512_cvtepi64_epi32(curve__compact1by2_avx512_64(_mm512_srli_epi64(vd, 1))));
        _mm256_storeu_si256((__m256i *)(void *)(z + i),
            _mm512_cvtepi64_epi32(curve__compact1by2_avx512_64(_mm512_srli_epi64(vd, 2))));
    }
#elif defined(CURVE_HAVE_AVX2)
    for (; i < n - (n & 3u); i += 4) {
        __m256i vd = _mm256_loadu_si256((const __m256i *)(const void *)(d + i));
        _mm_storeu_si128((__m128i *)(void *)(x + i),
            curve__pack64lo32_avx2(curve__compact1by2_avx2_64(vd)));
        _mm_storeu_si128((__m128i *)(void *)(y + i),
            curve__pack64lo32_avx2(curve__compact1by2_avx2_64(_mm256_srli_epi64(vd, 1))));
        _mm_storeu_si128((__m128i *)(void *)(z + i),
            curve__pack64lo32_avx2(curve__compact1by2_avx2_64(_mm256_srli_epi64(vd, 2))));
    }
#endif
    for (; i < n; ++i) curve_morton3_decode(d[i], x + i, y + i, z + i);
}

CURVE_INLINE void curve_hilbert_xy2d_batch(const uint32_t *CURVE_RESTRICT x,
                                           const uint32_t *CURVE_RESTRICT y,
                                           int bits, uint32_t *CURVE_RESTRICT out,
                                           size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) out[i] = curve_hilbert_xy2d(x[i], y[i], bits);
}

CURVE_INLINE void curve_hilbert_d2xy_batch(const uint32_t *CURVE_RESTRICT d,
                                           int bits,
                                           uint32_t *CURVE_RESTRICT x,
                                           uint32_t *CURVE_RESTRICT y,
                                           size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) curve_hilbert_d2xy(d[i], bits, x + i, y + i);
}

CURVE_INLINE void curve_hilbert2_64_xy2d_batch(const uint32_t *CURVE_RESTRICT x,
                                               const uint32_t *CURVE_RESTRICT y,
                                               int bits,
                                               uint64_t *CURVE_RESTRICT out,
                                               size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) out[i] = curve_hilbert2_64_xy2d(x[i], y[i], bits);
}

CURVE_INLINE void curve_hilbert2_64_d2xy_batch(const uint64_t *CURVE_RESTRICT d,
                                               int bits,
                                               uint32_t *CURVE_RESTRICT x,
                                               uint32_t *CURVE_RESTRICT y,
                                               size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) curve_hilbert2_64_d2xy(d[i], bits, x + i, y + i);
}

#ifdef __cplusplus
}
#endif

#endif /* CURVES_H */
