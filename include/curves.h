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
 *   • Wide n-D Morton / Hilbert keys (checked, caller-owned words, up to 256D)
 *   • Key prefix ops and spatial cell hierarchy (parent/child, AABB bounds)
 *   • AABB → key-interval decomposition (Morton / Hilbert) and prefix directory
 *   • Linear quantizers (float/double → bins) and fused quantize-to-key
 *   • Moore curve (cyclic Hilbert variant, bits 1..16)
 *   • Peano 2D (base-3 geometric curve, digits 1..10, invertible)
 *   • Gray-coded Morton 2D / 3D / n-D (Gray each axis, then interleave)
 *   • Mixed-bit (anisotropic) Morton wide keys
 *   • Sierpiński curve ordering key (Bartholdi & Platzman; non-dense, xy2d only)
 *   • Gilbert 2D / 3D visit (arbitrary rectangles, not just 2^n)
 *   • Batch APIs (restrict + NEON / AVX2 / AVX-512) for Gray, Morton, and
 *     Hilbert 2D (32-bit; 64-bit Hilbert on AVX2/AVX-512); compile-time
 *     ISA name via curve_compiled_isa / curve_simd_isa
 *   • Caller-owned runtime ISA dispatch (curve_dispatch), no process globals
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
 *   Define CURVE_RUNTIME_DISPATCH to compile extra x86 AVX2/AVX-512 batch
 *   kernels with target attributes for fat binaries (gcc/clang).
 *
 * Version: 1.1.0
 * License: MIT
 * ============================================================================ */

#ifndef CURVES_H
#define CURVES_H

#include <stdint.h>
#include <stddef.h>

#define CURVES_VERSION_MAJOR 1
#define CURVES_VERSION_MINOR 1
#define CURVES_VERSION_PATCH 0
#define CURVES_VERSION       "1.1.0"

#ifndef CURVE_MAX_DIMS
#define CURVE_MAX_DIMS 256
#endif
#define CURVE_MAX_BITS_PER_AXIS 32u

/* Box → interval decomposition enumerates 2^n children per cell. */
#ifndef CURVE_RANGE_MAX_DIMS
#define CURVE_RANGE_MAX_DIMS 8
#endif

/* -------------------------------------------------------------------------- */
/* Compiler / ISA helpers                                                     */
/* -------------------------------------------------------------------------- */

#if defined(CURVE_RUNTIME_DISPATCH) && !defined(CURVE_NO_SIMD) && \
    (defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)) && \
    (defined(__GNUC__) || defined(__clang__))
#  define CURVE_DISPATCH_X86 1
#endif

#if defined(__AVX512F__) || defined(__AVX2__) || defined(__BMI2__) || defined(CURVE_DISPATCH_X86)
#  include <immintrin.h>
#endif
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
#  if defined(_MSC_VER)
#    include <intrin.h>
#  elif defined(__GNUC__) || defined(__clang__)
#    include <cpuid.h>
#  endif
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

#if defined(CURVE_DISPATCH_X86) && !defined(CURVE_HAVE_AVX2)
#  define CURVE_ATTR_AVX2 __attribute__((target("avx2")))
#else
#  define CURVE_ATTR_AVX2
#endif
#if defined(CURVE_DISPATCH_X86) && !defined(CURVE_HAVE_AVX512)
#  define CURVE_ATTR_AVX512 __attribute__((target("avx512f")))
#else
#  define CURVE_ATTR_AVX512
#endif

#define CURVE_ISA_SCALAR  (1u << 0)
#define CURVE_ISA_NEON    (1u << 1)
#define CURVE_ISA_AVX2    (1u << 2)
#define CURVE_ISA_AVX512  (1u << 3)
#define CURVE_ISA_BMI2    (1u << 4)

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
#  define CURVE__SA_JOIN(a, b) a##b
#  define CURVE__SA_X(a, b)    CURVE__SA_JOIN(a, b)
#  define CURVE_STATIC_ASSERT(cond, msg) \
        struct CURVE__SA_X(curve__sa_, __LINE__) { char curve__sa_c; }
#endif

/* Checked-API status. Legacy encode/decode functions keep their current
 * return types and still clamp rather than report errors. */
typedef enum curve_status {
    CURVE_OK           = 0,
    CURVE_EINVAL       = 1, /* null, bad dims/bits/bounds, noncanonical, partial cell */
    CURVE_ERANGE       = 2, /* dims > MAX, no successor/parent, leaf has no child */
    CURVE_ENOSPACE     = 3, /* output or scratch capacity too small */
    CURVE_EOVERFLOW    = 4, /* size arithmetic overflow */
    CURVE_EUNSUPPORTED = 5  /* child index n_dims>64, box walk n_dims too large */
} curve_status;

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

/* -------------------------------------------------------------------------- */
/* Gray-coded Morton — binary-reflected Gray on each axis, then interleave    */
/*                                                                            */
/* Domain matches the corresponding Morton map (16 / 32 / 21 bits per axis).  */
/* Mask first, then Gray-encode; inverse is compact then Gray-decode.         */
/* -------------------------------------------------------------------------- */

CURVE_INLINE uint32_t curve__gray32(uint32_t x)
{
    return x ^ (x >> 1);
}

CURVE_INLINE uint32_t curve__ungray32(uint32_t g)
{
    return (uint32_t)curve_gray_decode((uint64_t)g);
}

CURVE_INLINE uint32_t curve_gmorton2(uint32_t x, uint32_t y)
{
    return curve_morton2(curve__gray32(x & 0xFFFFu), curve__gray32(y & 0xFFFFu));
}

CURVE_INLINE void curve_gmorton2_decode(uint32_t d, uint32_t *CURVE_RESTRICT x,
                                        uint32_t *CURVE_RESTRICT y)
{
    uint32_t gx, gy;
    curve_morton2_decode(d, &gx, &gy);
    *x = curve__ungray32(gx);
    *y = curve__ungray32(gy);
}

CURVE_INLINE uint64_t curve_gmorton2_64(uint32_t x, uint32_t y)
{
    return curve_morton2_64(curve__gray32(x), curve__gray32(y));
}

CURVE_INLINE void curve_gmorton2_64_decode(uint64_t d, uint32_t *CURVE_RESTRICT x,
                                           uint32_t *CURVE_RESTRICT y)
{
    uint32_t gx, gy;
    curve_morton2_64_decode(d, &gx, &gy);
    *x = curve__ungray32(gx);
    *y = curve__ungray32(gy);
}

CURVE_INLINE uint64_t curve_gmorton3(uint32_t x, uint32_t y, uint32_t z)
{
    return curve_morton3(curve__gray32(x & 0x1FFFFFu),
                         curve__gray32(y & 0x1FFFFFu),
                         curve__gray32(z & 0x1FFFFFu));
}

CURVE_INLINE void curve_gmorton3_decode(uint64_t d, uint32_t *CURVE_RESTRICT x,
                                        uint32_t *CURVE_RESTRICT y,
                                        uint32_t *CURVE_RESTRICT z)
{
    uint32_t gx, gy, gz;
    curve_morton3_decode(d, &gx, &gy, &gz);
    *x = curve__ungray32(gx);
    *y = curve__ungray32(gy);
    *z = curve__ungray32(gz);
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
/* Peano 2D — geometric base-3 curve, digits in 1..10                         */
/*                                                                            */
/* Grid 3^digits × 3^digits, index 0 .. 9^digits - 1 (fits uint32).           */
/* Generator (origin bottom-left, x right, y up) is the column snake:         */
/*   y=2:  2 3 8                                                              */
/*   y=1:  1 4 7                                                              */
/*   y=0:  0 5 6                                                              */
/* Sub-squares of order n-1 are reflected (Rx / Ry / R180) so consecutive     */
/* indices are edge-adjacent.  This is the connected geometric Peano curve,   */
/* not ternary-digit interleave (Lebesgue in base 3).                         */
/* -------------------------------------------------------------------------- */

CURVE_ALIGNAS(64) static const uint32_t curve__pow3[11] = {
    1u, 3u, 9u, 27u, 81u, 243u, 729u, 2187u, 6561u, 19683u, 59049u
};
static const uint8_t curve__peano_q2x[9] = { 0, 0, 0, 1, 1, 1, 2, 2, 2 };
static const uint8_t curve__peano_q2y[9] = { 0, 1, 2, 2, 1, 0, 0, 1, 2 };
static const uint8_t curve__peano_xy2q[3][3] = {
    { 0, 1, 2 },
    { 5, 4, 3 },
    { 6, 7, 8 }
};
/* Per-cell remainder reflection: bit0 = flip x, bit1 = flip y. */
static const uint8_t curve__peano_tf[9] = { 0, 1, 0, 2, 3, 2, 0, 1, 0 };

CURVE_INLINE int curve__clamp_peano_digits(int digits)
{
    if (digits < 1)  return 1;
    if (digits > 10) return 10;
    return digits;
}

CURVE_INLINE uint32_t curve_peano_xy2d(uint32_t x, uint32_t y, int digits)
{
    uint32_t d = 0, p, n;
    digits = curve__clamp_peano_digits(digits);
    n = curve__pow3[digits];
    x %= n;
    y %= n;
    for (p = curve__pow3[digits - 1]; p > 0u; p /= 3u) {
        uint32_t xi = x / p, yi = y / p, q, tf;
        x %= p;
        y %= p;
        q = curve__peano_xy2q[xi][yi];
        d = d * 9u + q;
        tf = curve__peano_tf[q];
        if (tf & 1u) x = (p - 1u) - x;
        if (tf & 2u) y = (p - 1u) - y;
    }
    return d;
}

CURVE_INLINE void curve_peano_d2xy(uint32_t d, int digits,
                                   uint32_t *CURVE_RESTRICT x,
                                   uint32_t *CURVE_RESTRICT y)
{
    uint32_t xx = 0, yy = 0, div = 1u;
    unsigned t = 0;
    int k;
    digits = curve__clamp_peano_digits(digits);
    for (k = 1; k < digits; k++) div *= 9u;
    d %= curve__pow3[digits] * curve__pow3[digits];
    for (k = 0; k < digits; k++) {
        uint32_t q = d / div;
        uint32_t lx = curve__peano_q2x[q], ly = curve__peano_q2y[q];
        uint32_t gx = (t & 1u) ? 2u - lx : lx;
        uint32_t gy = (t & 2u) ? 2u - ly : ly;
        d %= div;
        div /= 9u;
        xx = xx * 3u + gx;
        yy = yy * 3u + gy;
        t ^= curve__peano_tf[q];
    }
    *x = xx;
    *y = yy;
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
/* Wide n-D keys — canonical words / bytes, checked Morton and Hilbert        */
/*                                                                            */
/* A key has K = dimensions * bits_per_axis logical bits.                     */
/*                                                                            */
/* Word view: W = ceil(K/64) uint64_t words; words[0] holds the most          */
/* significant 64 bits. Unused bits of the last word are left-aligned         */
/* (significant bits in the high end, unused low bits zero).                  */
/* When K <= 64 the legacy integer index h satisfies                          */
/*   words[0] == h << (64-K).                                                 */
/*                                                                            */
/* Byte view: B = ceil(K/8) bytes, big-endian of the logical bits; unused     */
/* low bits of the last byte are zero. Word/byte conversion is independent    */
/* of host endianness.                                                        */
/*                                                                            */
/* Packing (Morton and Hilbert): for bit plane p = bits-1 .. 0, emit axis     */
/* 0 then 1 then … (axis 0 is more significant within a level). This matches  */
/* curve_hilbertnd_xy2d. It does NOT match curve_morton2 / curve_morton3,     */
/* which put axis 0 on even bits counting from the LSB.                       */
/*                                                                            */
/* Checked APIs validate every argument before writing outputs.               */
/* dimensions == 0 or bits outside 1..32 → CURVE_EINVAL (no clamp).           */
/* Hilbert requires a caller-owned scratch buffer of `dimensions` uint32_t;   */
/* scratch may alias coords (in-place) but must not alias key_words.          */
/* Morton does not use scratch.                                               */
/* -------------------------------------------------------------------------- */

CURVE_INLINE curve_status curve__check_dims_bits(size_t dimensions,
                                                 unsigned bits_per_axis)
{
    if (dimensions < 1u || bits_per_axis < 1u ||
        bits_per_axis > CURVE_MAX_BITS_PER_AXIS)
        return CURVE_EINVAL;
    if (dimensions > (size_t)CURVE_MAX_DIMS)
        return CURVE_ERANGE;
    return CURVE_OK;
}

CURVE_INLINE size_t curve_key_bits(size_t dimensions, unsigned bits_per_axis)
{
    if (curve__check_dims_bits(dimensions, bits_per_axis) != CURVE_OK)
        return 0;
    if (bits_per_axis != 0u && dimensions > ((size_t)-1) / bits_per_axis)
        return 0;
    return dimensions * (size_t)bits_per_axis;
}

CURVE_INLINE size_t curve_key_words_required(size_t dimensions,
                                             unsigned bits_per_axis)
{
    size_t bits = curve_key_bits(dimensions, bits_per_axis);
    if (bits == 0) return 0;
    return (bits + 63u) / 64u;
}

CURVE_INLINE size_t curve_key_bytes_required(size_t dimensions,
                                             unsigned bits_per_axis)
{
    size_t bits = curve_key_bits(dimensions, bits_per_axis);
    if (bits == 0) return 0;
    return (bits + 7u) / 8u;
}

CURVE_INLINE unsigned curve__key_last_word_unused(size_t key_bits)
{
    return (unsigned)((64u - (key_bits % 64u)) % 64u);
}

CURVE_INLINE uint64_t curve__low_ones64(unsigned n)
{
    if (n == 0u) return 0ull;
    if (n >= 64u) return ~0ull;
    return (1ull << n) - 1ull;
}

CURVE_INLINE curve_status curve__key_check_canonical(const uint64_t *words,
                                                     size_t n_words,
                                                     size_t key_bits)
{
    unsigned unused;
    size_t need = (key_bits + 63u) / 64u;
    if (key_bits == 0u || words == NULL || n_words != need)
        return CURVE_EINVAL;
    unused = curve__key_last_word_unused(key_bits);
    if (unused && (words[need - 1u] & curve__low_ones64(unused)) != 0ull)
        return CURVE_EINVAL;
    return CURVE_OK;
}

CURVE_INLINE curve_status curve_key_words_to_bytes(const uint64_t *words,
                                                   size_t key_bits,
                                                   uint8_t *bytes,
                                                   size_t byte_capacity)
{
    size_t n_bytes, i;
    unsigned unused_byte_bits;
    if (key_bits == 0u || words == NULL || bytes == NULL)
        return CURVE_EINVAL;
    n_bytes = (key_bits + 7u) / 8u;
    if (byte_capacity < n_bytes) return CURVE_ENOSPACE;
    for (i = 0; i < n_bytes; i++) {
        size_t bit = i * 8u;
        size_t wi = bit / 64u;
        unsigned shift = (unsigned)(56u - (bit % 64u));
        bytes[i] = (uint8_t)(words[wi] >> shift);
    }
    unused_byte_bits = (unsigned)((8u - (key_bits % 8u)) % 8u);
    if (unused_byte_bits)
        bytes[n_bytes - 1u] = (uint8_t)(bytes[n_bytes - 1u] &
            (uint8_t)(0xFFu << unused_byte_bits));
    return CURVE_OK;
}

CURVE_INLINE curve_status curve_key_bytes_to_words(const uint8_t *bytes,
                                                   size_t key_bits,
                                                   uint64_t *words,
                                                   size_t word_capacity)
{
    size_t n_bytes, n_words, i;
    unsigned unused_byte_bits, unused_word_bits;
    if (key_bits == 0u || bytes == NULL || words == NULL)
        return CURVE_EINVAL;
    n_bytes = (key_bits + 7u) / 8u;
    n_words = (key_bits + 63u) / 64u;
    if (word_capacity < n_words) return CURVE_ENOSPACE;
    unused_byte_bits = (unsigned)((8u - (key_bits % 8u)) % 8u);
    if (unused_byte_bits &&
        (bytes[n_bytes - 1u] & (uint8_t)((1u << unused_byte_bits) - 1u)) != 0u)
        return CURVE_EINVAL;
    for (i = 0; i < n_words; i++) words[i] = 0;
    for (i = 0; i < n_bytes; i++) {
        size_t bit = i * 8u;
        size_t wi = bit / 64u;
        unsigned shift = (unsigned)(56u - (bit % 64u));
        words[wi] |= ((uint64_t)bytes[i]) << shift;
    }
    unused_word_bits = curve__key_last_word_unused(key_bits);
    if (unused_word_bits)
        words[n_words - 1u] &= ~curve__low_ones64(unused_word_bits);
    return CURVE_OK;
}

CURVE_INLINE int curve_key_compare(const uint64_t *a, const uint64_t *b,
                                   size_t n_words)
{
    size_t i;
    if (n_words == 0u) return 0;
    for (i = 0; i < n_words; i++) {
        if (a[i] < b[i]) return -1;
        if (a[i] > b[i]) return 1;
    }
    return 0;
}

CURVE_INLINE int curve_key_equal(const uint64_t *a, const uint64_t *b,
                                 size_t n_words)
{
    return curve_key_compare(a, b, n_words) == 0;
}

/* -------------------------------------------------------------------------- */
/* Prefix operations on canonical wide keys                                   */
/*                                                                            */
/* Bit 0 is the most significant logical bit of the key (MSB of words[0]).    */
/* A prefix of length P is those first P bits. Remaining bits of a masked     */
/* prefix are zero. The prefix interval is [masked, successor) when a         */
/* successor exists; otherwise [masked, +∞) with out_unbounded set.           */
/* prefix_bits == 0 covers the whole domain and has no successor.             */
/* An all-ones prefix has no successor. CURVE_ERANGE on successor does not    */
/* write `out`. No 64-bit shift of a 64-bit value.                            */
/* -------------------------------------------------------------------------- */

CURVE_INLINE unsigned curve__clz64(uint64_t x)
{
#if defined(__GNUC__) || defined(__clang__)
    return (unsigned)__builtin_clzll(x);
#else
    unsigned n = 0;
    if ((x >> 32) == 0ull) { n += 32u; x <<= 32; }
    if ((x >> 48) == 0ull) { n += 16u; x <<= 16; }
    if ((x >> 56) == 0ull) { n += 8u;  x <<= 8; }
    if ((x >> 60) == 0ull) { n += 4u;  x <<= 4; }
    if ((x >> 62) == 0ull) { n += 2u;  x <<= 2; }
    if ((x >> 63) == 0ull) { n += 1u; }
    return n;
#endif
}

CURVE_INLINE size_t curve__key_n_words(size_t key_bits)
{
    return (key_bits + 63u) / 64u;
}

/* Bit 0 is the most significant logical bit of the key. */
CURVE_INLINE unsigned curve__key_get_bit(const uint64_t *words,
                                         size_t bit_from_msb)
{
    size_t wi = bit_from_msb / 64u;
    unsigned sh = (unsigned)(63u - (bit_from_msb % 64u));
    return (unsigned)((words[wi] >> sh) & 1ull);
}

CURVE_INLINE void curve__key_set_bit(uint64_t *words, size_t bit_from_msb,
                                     unsigned bit)
{
    size_t wi = bit_from_msb / 64u;
    unsigned sh = (unsigned)(63u - (bit_from_msb % 64u));
    uint64_t m = 1ull << sh;
    if (bit)
        words[wi] |= m;
    else
        words[wi] &= ~m;
}

CURVE_INLINE curve_status curve_key_get_bit(const uint64_t *words,
                                            size_t key_bits,
                                            size_t bit_from_msb,
                                            unsigned *out_bit)
{
    if (words == NULL || out_bit == NULL || key_bits == 0u)
        return CURVE_EINVAL;
    if (bit_from_msb >= key_bits) return CURVE_EINVAL;
    *out_bit = curve__key_get_bit(words, bit_from_msb);
    return CURVE_OK;
}

CURVE_INLINE curve_status curve_key_common_prefix_bits(const uint64_t *a,
                                                       const uint64_t *b,
                                                       size_t key_bits,
                                                       size_t *out_prefix)
{
    size_t W, i, used;
    if (a == NULL || b == NULL || out_prefix == NULL || key_bits == 0u)
        return CURVE_EINVAL;
    W = curve__key_n_words(key_bits);
    for (i = 0; i < W; i++) {
        uint64_t x = a[i] ^ b[i];
        used = ((i + 1u) * 64u <= key_bits) ? 64u : (key_bits - i * 64u);
        if (used < 64u)
            x &= ~curve__low_ones64((unsigned)(64u - used));
        if (x != 0ull) {
            unsigned lz = curve__clz64(x);
            size_t p = i * 64u + (size_t)lz;
            if (p > key_bits) p = key_bits;
            *out_prefix = p;
            return CURVE_OK;
        }
    }
    *out_prefix = key_bits;
    return CURVE_OK;
}

CURVE_INLINE curve_status curve_key_prefix_mask(const uint64_t *in,
                                                size_t key_bits,
                                                size_t prefix_bits,
                                                uint64_t *out, size_t out_words)
{
    size_t W, i, full, wi;
    unsigned keep;
    if (in == NULL || out == NULL || key_bits == 0u)
        return CURVE_EINVAL;
    if (prefix_bits > key_bits) return CURVE_EINVAL;
    W = curve__key_n_words(key_bits);
    if (out_words < W) return CURVE_ENOSPACE;
    if (prefix_bits == 0u) {
        for (i = 0; i < W; i++) out[i] = 0;
        return CURVE_OK;
    }
    for (i = 0; i < W; i++) out[i] = in[i];
    full = prefix_bits / 64u;
    keep = (unsigned)(prefix_bits % 64u);
    if (keep == 0u) {
        for (i = full; i < W; i++) out[i] = 0;
        return CURVE_OK;
    }
    wi = full;
    out[wi] &= ~curve__low_ones64((unsigned)(64u - keep));
    for (i = wi + 1u; i < W; i++) out[i] = 0;
    return CURVE_OK;
}

CURVE_INLINE curve_status curve_key_prefix_matches(const uint64_t *key,
                                                   const uint64_t *prefix,
                                                   size_t key_bits,
                                                   size_t prefix_bits,
                                                   int *out_match)
{
    size_t i, full;
    unsigned keep;
    uint64_t mask;
    if (key == NULL || prefix == NULL || out_match == NULL || key_bits == 0u)
        return CURVE_EINVAL;
    if (prefix_bits > key_bits) return CURVE_EINVAL;
    if (prefix_bits == 0u) {
        *out_match = 1;
        return CURVE_OK;
    }
    full = prefix_bits / 64u;
    for (i = 0; i < full; i++) {
        if (key[i] != prefix[i]) {
            *out_match = 0;
            return CURVE_OK;
        }
    }
    keep = (unsigned)(prefix_bits % 64u);
    if (keep == 0u) {
        *out_match = 1;
        return CURVE_OK;
    }
    mask = ~curve__low_ones64((unsigned)(64u - keep));
    *out_match = ((key[full] & mask) == (prefix[full] & mask));
    return CURVE_OK;
}

/* Increment the prefix_bits high bits by one. Remaining bits of `out` are
 * zero. CURVE_ERANGE means the prefix is all ones (no successor): the
 * interval is [masked prefix, +∞). prefix_bits == 0 also has no successor
 * (the prefix covers the whole domain). */
CURVE_INLINE curve_status curve_key_prefix_successor(const uint64_t *in,
                                                     size_t key_bits,
                                                     size_t prefix_bits,
                                                     uint64_t *out,
                                                     size_t out_words)
{
    curve_status st;
    size_t bit, i, full;
    unsigned keep;
    int carry, all_ones;
    if (prefix_bits == 0u) return CURVE_ERANGE;
    if (in == NULL || key_bits == 0u) return CURVE_EINVAL;
    if (prefix_bits > key_bits) return CURVE_EINVAL;
    full = prefix_bits / 64u;
    keep = (unsigned)(prefix_bits % 64u);
    all_ones = 1;
    for (i = 0; i < full; i++) {
        if (in[i] != ~0ull) {
            all_ones = 0;
            break;
        }
    }
    if (all_ones && keep == 0u) return CURVE_ERANGE;
    if (all_ones && keep != 0u) {
        uint64_t hi = ~curve__low_ones64((unsigned)(64u - keep));
        if ((in[full] & hi) == hi) return CURVE_ERANGE;
    }
    st = curve_key_prefix_mask(in, key_bits, prefix_bits, out, out_words);
    if (st != CURVE_OK) return st;
    carry = 1;
    bit = prefix_bits;
    while (bit > 0u && carry) {
        unsigned cur;
        bit--;
        cur = curve__key_get_bit(out, bit);
        curve__key_set_bit(out, bit, cur ? 0u : 1u);
        carry = cur ? 1 : 0;
    }
    return CURVE_OK;
}

CURVE_INLINE curve_status curve_key_prefix_interval(const uint64_t *in,
                                                    size_t key_bits,
                                                    size_t prefix_bits,
                                                    uint64_t *lower,
                                                    uint64_t *upper,
                                                    size_t out_words,
                                                    int *out_unbounded)
{
    curve_status st;
    if (lower == NULL || upper == NULL || out_unbounded == NULL)
        return CURVE_EINVAL;
    st = curve_key_prefix_mask(in, key_bits, prefix_bits, lower, out_words);
    if (st != CURVE_OK) return st;
    st = curve_key_prefix_successor(in, key_bits, prefix_bits, upper,
                                    out_words);
    if (st == CURVE_ERANGE) {
        size_t i, W = curve__key_n_words(key_bits);
        if (out_words < W) return CURVE_ENOSPACE;
        for (i = 0; i < W; i++) upper[i] = 0;
        *out_unbounded = 1;
        return CURVE_OK;
    }
    if (st != CURVE_OK) return st;
    *out_unbounded = 0;
    return CURVE_OK;
}

CURVE_INLINE size_t curve_key_cell_level(size_t prefix_bits, size_t n_dims)
{
    if (n_dims == 0u) return 0;
    return prefix_bits / n_dims;
}

CURVE_INLINE size_t curve_key_cell_prefix_bits(size_t level, size_t n_dims)
{
    if (n_dims != 0u && level > ((size_t)-1) / n_dims) return 0;
    return level * n_dims;
}

CURVE_INLINE int curve_key_prefix_is_cell(size_t prefix_bits, size_t n_dims,
                                          unsigned bits_per_axis)
{
    size_t max_bits;
    if (n_dims == 0u || bits_per_axis < 1u ||
        bits_per_axis > CURVE_MAX_BITS_PER_AXIS)
        return 0;
    if (prefix_bits % n_dims != 0u) return 0;
    max_bits = n_dims * (size_t)bits_per_axis;
    return prefix_bits <= max_bits;
}

/* First index i in a sorted array of n_keys keys (stride `stride_words`)
 * such that keys[i] >= query. Returns n_keys if every key is < query.
 * n_words is the comparison width; stride_words must be >= n_words. */
CURVE_INLINE curve_status curve_key_lower_bound(const uint64_t *keys,
                                                size_t n_keys,
                                                size_t stride_words,
                                                const uint64_t *query,
                                                size_t n_words,
                                                size_t *out_index)
{
    size_t lo, hi;
    if (out_index == NULL || n_words == 0u) return CURVE_EINVAL;
    if (n_keys != 0u && (keys == NULL || query == NULL)) return CURVE_EINVAL;
    if (n_keys != 0u && stride_words < n_words) return CURVE_EINVAL;
    if (n_keys == 0u) {
        *out_index = 0;
        return CURVE_OK;
    }
    lo = 0;
    hi = n_keys;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2u;
        const uint64_t *k = keys + mid * stride_words;
        if (n_words == 1u) {
            if (k[0] < query[0]) lo = mid + 1u;
            else hi = mid;
        } else if (curve_key_compare(k, query, n_words) < 0) {
            lo = mid + 1u;
        } else {
            hi = mid;
        }
    }
    *out_index = lo;
    return CURVE_OK;
}

CURVE_INLINE curve_status curve_key_upper_bound(const uint64_t *keys,
                                                size_t n_keys,
                                                size_t stride_words,
                                                const uint64_t *query,
                                                size_t n_words,
                                                size_t *out_index)
{
    size_t lo, hi;
    if (out_index == NULL || n_words == 0u) return CURVE_EINVAL;
    if (n_keys != 0u && (keys == NULL || query == NULL)) return CURVE_EINVAL;
    if (n_keys != 0u && stride_words < n_words) return CURVE_EINVAL;
    if (n_keys == 0u) {
        *out_index = 0;
        return CURVE_OK;
    }
    lo = 0;
    hi = n_keys;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2u;
        const uint64_t *k = keys + mid * stride_words;
        int cmp;
        if (n_words == 1u)
            cmp = (k[0] < query[0]) ? -1 : (k[0] > query[0] ? 1 : 0);
        else
            cmp = curve_key_compare(k, query, n_words);
        if (cmp <= 0) lo = mid + 1u;
        else hi = mid;
    }
    *out_index = lo;
    return CURVE_OK;
}

/* [begin, end) of keys sharing `prefix_bits` of `prefix`. end == n_keys
 * when the prefix has no successor. `scratch` must hold 2 * n_words
 * uint64_t (lower then upper interval endpoints). */
CURVE_INLINE curve_status curve_key_prefix_range(const uint64_t *keys,
                                                 size_t n_keys,
                                                 size_t stride_words,
                                                 const uint64_t *prefix,
                                                 size_t key_bits,
                                                 size_t prefix_bits,
                                                 uint64_t *scratch,
                                                 size_t scratch_words,
                                                 size_t *out_begin,
                                                 size_t *out_end)
{
    curve_status st;
    size_t W, begin, end;
    uint64_t *lower, *upper;
    int unbounded = 0;
    if (out_begin == NULL || out_end == NULL || scratch == NULL)
        return CURVE_EINVAL;
    if (key_bits == 0u) return CURVE_EINVAL;
    W = curve__key_n_words(key_bits);
    if (scratch_words < 2u * W) return CURVE_ENOSPACE;
    lower = scratch;
    upper = scratch + W;
    st = curve_key_prefix_interval(prefix, key_bits, prefix_bits,
                                   lower, upper, W, &unbounded);
    if (st != CURVE_OK) return st;
    st = curve_key_lower_bound(keys, n_keys, stride_words, lower, W, &begin);
    if (st != CURVE_OK) return st;
    if (unbounded)
        end = n_keys;
    else {
        st = curve_key_lower_bound(keys, n_keys, stride_words, upper, W, &end);
        if (st != CURVE_OK) return st;
    }
    *out_begin = begin;
    *out_end = end;
    return CURVE_OK;
}

/* In-place Skilling transform. n_dims == 1 is a no-op. Coordinates are not
 * masked; callers that need a `bits`-bit grid should mask first. */
CURVE_INLINE curve_status curve_axes_to_transpose(uint32_t *axes, size_t n_dims,
                                                  unsigned bits)
{
    curve_status st = curve__check_dims_bits(n_dims, bits);
    if (st != CURVE_OK) return st;
    if (axes == NULL) return CURVE_EINVAL;
    if (n_dims == 1u) return CURVE_OK;
    curve__axes_to_transpose(axes, (int)bits, (int)n_dims);
    return CURVE_OK;
}

CURVE_INLINE curve_status curve_transpose_to_axes(uint32_t *axes, size_t n_dims,
                                                  unsigned bits)
{
    curve_status st = curve__check_dims_bits(n_dims, bits);
    if (st != CURVE_OK) return st;
    if (axes == NULL) return CURVE_EINVAL;
    if (n_dims == 1u) return CURVE_OK;
    curve__transpose_to_axes(axes, (int)bits, (int)n_dims);
    return CURVE_OK;
}

CURVE_INLINE curve_status curve__pack_axes_words(const uint32_t *axes, size_t n,
                                                 unsigned bits, uint64_t *words,
                                                 size_t cap)
{
    size_t K, W, k, i;
    unsigned p;
    uint32_t mask;
    K = n * (size_t)bits;
    W = (K + 63u) / 64u;
    if (cap < W) return CURVE_ENOSPACE;
    mask = curve__lo_mask((int)bits);
    for (i = 0; i < W; i++) words[i] = 0;
    k = 0;
    p = bits;
    while (p-- > 0u) {
        for (i = 0; i < n; i++) {
            uint64_t bit = ((uint64_t)(axes[i] & mask) >> p) & 1ull;
            words[k / 64u] |= bit << (63u - (unsigned)(k % 64u));
            k++;
        }
    }
    return CURVE_OK;
}

CURVE_INLINE curve_status curve__unpack_axes_words(const uint64_t *words,
                                                   size_t n_words, size_t n,
                                                   unsigned bits,
                                                   uint32_t *axes, size_t cap)
{
    size_t K, k, i;
    unsigned p;
    curve_status st;
    K = n * (size_t)bits;
    st = curve__key_check_canonical(words, n_words, K);
    if (st != CURVE_OK) return st;
    if (axes == NULL) return CURVE_EINVAL;
    if (cap < n) return CURVE_ENOSPACE;
    for (i = 0; i < n; i++) axes[i] = 0;
    k = 0;
    p = bits;
    while (p-- > 0u) {
        for (i = 0; i < n; i++) {
            uint32_t bit = (uint32_t)((words[k / 64u] >>
                (63u - (unsigned)(k % 64u))) & 1ull);
            axes[i] |= bit << p;
            k++;
        }
    }
    return CURVE_OK;
}

CURVE_INLINE curve_status curve_mortonnd_encode_words(
    const uint32_t *coords, size_t dimensions, unsigned bits_per_axis,
    uint64_t *key_words, size_t key_capacity_words)
{
    curve_status st = curve__check_dims_bits(dimensions, bits_per_axis);
    size_t need;
    if (st != CURVE_OK) return st;
    if (coords == NULL || key_words == NULL) return CURVE_EINVAL;
    need = curve_key_words_required(dimensions, bits_per_axis);
    if (need == 0u) return CURVE_EOVERFLOW;
    if (key_capacity_words < need) return CURVE_ENOSPACE;
    return curve__pack_axes_words(coords, dimensions, bits_per_axis,
                                  key_words, key_capacity_words);
}

CURVE_INLINE curve_status curve_mortonnd_decode_words(
    const uint64_t *key_words, size_t key_count_words,
    size_t dimensions, unsigned bits_per_axis,
    uint32_t *coords, size_t coords_capacity)
{
    curve_status st = curve__check_dims_bits(dimensions, bits_per_axis);
    if (st != CURVE_OK) return st;
    if (key_words == NULL || coords == NULL) return CURVE_EINVAL;
    if (curve_key_words_required(dimensions, bits_per_axis) == 0u)
        return CURVE_EOVERFLOW;
    return curve__unpack_axes_words(key_words, key_count_words, dimensions,
                                    bits_per_axis, coords, coords_capacity);
}

CURVE_INLINE curve_status curve_hilbertnd_encode_words(
    const uint32_t *coords, size_t dimensions, unsigned bits_per_axis,
    uint64_t *key_words, size_t key_capacity_words,
    uint32_t *scratch_axes, size_t scratch_capacity_axes)
{
    curve_status st = curve__check_dims_bits(dimensions, bits_per_axis);
    size_t need, i;
    uint32_t mask;
    if (st != CURVE_OK) return st;
    if (coords == NULL || key_words == NULL || scratch_axes == NULL)
        return CURVE_EINVAL;
    need = curve_key_words_required(dimensions, bits_per_axis);
    if (need == 0u) return CURVE_EOVERFLOW;
    if (key_capacity_words < need) return CURVE_ENOSPACE;
    if (scratch_capacity_axes < dimensions) return CURVE_ENOSPACE;
    mask = curve__lo_mask((int)bits_per_axis);
    for (i = 0; i < dimensions; i++)
        scratch_axes[i] = coords[i] & mask;
    st = curve_axes_to_transpose(scratch_axes, dimensions, bits_per_axis);
    if (st != CURVE_OK) return st;
    return curve__pack_axes_words(scratch_axes, dimensions, bits_per_axis,
                                  key_words, key_capacity_words);
}

CURVE_INLINE curve_status curve_hilbertnd_decode_words(
    const uint64_t *key_words, size_t key_count_words,
    size_t dimensions, unsigned bits_per_axis,
    uint32_t *coords, size_t coords_capacity,
    uint32_t *scratch_axes, size_t scratch_capacity_axes)
{
    curve_status st = curve__check_dims_bits(dimensions, bits_per_axis);
    size_t i;
    if (st != CURVE_OK) return st;
    if (key_words == NULL || coords == NULL || scratch_axes == NULL)
        return CURVE_EINVAL;
    if (coords_capacity < dimensions) return CURVE_ENOSPACE;
    if (scratch_capacity_axes < dimensions) return CURVE_ENOSPACE;
    st = curve__unpack_axes_words(key_words, key_count_words, dimensions,
                                  bits_per_axis, scratch_axes,
                                  scratch_capacity_axes);
    if (st != CURVE_OK) return st;
    st = curve_transpose_to_axes(scratch_axes, dimensions, bits_per_axis);
    if (st != CURVE_OK) return st;
    if (scratch_axes != coords) {
        for (i = 0; i < dimensions; i++)
            coords[i] = scratch_axes[i];
    }
    return CURVE_OK;
}

CURVE_INLINE curve_status curve_gmortonnd_encode_words(
    const uint32_t *coords, size_t dimensions, unsigned bits_per_axis,
    uint64_t *key_words, size_t key_capacity_words)
{
    curve_status st = curve__check_dims_bits(dimensions, bits_per_axis);
    size_t need, i, K, W, k;
    unsigned p;
    uint32_t mask;
    if (st != CURVE_OK) return st;
    if (coords == NULL || key_words == NULL) return CURVE_EINVAL;
    need = curve_key_words_required(dimensions, bits_per_axis);
    if (need == 0u) return CURVE_EOVERFLOW;
    if (key_capacity_words < need) return CURVE_ENOSPACE;
    mask = curve__lo_mask((int)bits_per_axis);
    K = dimensions * (size_t)bits_per_axis;
    W = (K + 63u) / 64u;
    for (i = 0; i < W; i++) key_words[i] = 0;
    k = 0;
    p = bits_per_axis;
    while (p-- > 0u) {
        for (i = 0; i < dimensions; i++) {
            uint32_t g = curve__gray32(coords[i] & mask);
            uint64_t bit = ((uint64_t)g >> p) & 1ull;
            key_words[k / 64u] |= bit << (63u - (unsigned)(k % 64u));
            k++;
        }
    }
    return CURVE_OK;
}

CURVE_INLINE curve_status curve_gmortonnd_decode_words(
    const uint64_t *key_words, size_t key_count_words,
    size_t dimensions, unsigned bits_per_axis,
    uint32_t *coords, size_t coords_capacity)
{
    curve_status st;
    size_t i;
    st = curve_mortonnd_decode_words(key_words, key_count_words, dimensions,
                                     bits_per_axis, coords, coords_capacity);
    if (st != CURVE_OK) return st;
    for (i = 0; i < dimensions; i++)
        coords[i] = curve__ungray32(coords[i]);
    return CURVE_OK;
}

/* Mixed-bit (anisotropic) Morton. Each axis has its own 1..32 bit width.
 * Planes are emitted from max(bits)-1 down to 0; within a plane, axis 0
 * first, skipping axes that do not have that plane. The result is a
 * rectangular product of dyadic intervals, not a Hilbert curve. Uniform
 * bits match curve_mortonnd_encode_words. */

CURVE_INLINE curve_status curve__check_mixed_bits(size_t dimensions,
                                                  const unsigned *bits_per_axis)
{
    size_t i;
    if (dimensions < 1u) return CURVE_EINVAL;
    if (dimensions > (size_t)CURVE_MAX_DIMS) return CURVE_ERANGE;
    if (bits_per_axis == NULL) return CURVE_EINVAL;
    for (i = 0; i < dimensions; i++) {
        if (bits_per_axis[i] < 1u ||
            bits_per_axis[i] > CURVE_MAX_BITS_PER_AXIS)
            return CURVE_EINVAL;
    }
    return CURVE_OK;
}

CURVE_INLINE size_t curve_key_bits_mixed(size_t dimensions,
                                         const unsigned *bits_per_axis)
{
    size_t i, s = 0;
    if (curve__check_mixed_bits(dimensions, bits_per_axis) != CURVE_OK)
        return 0;
    for (i = 0; i < dimensions; i++) {
        if (s > (size_t)-1 - (size_t)bits_per_axis[i])
            return 0;
        s += (size_t)bits_per_axis[i];
    }
    return s;
}

CURVE_INLINE size_t curve_key_words_required_mixed(
    size_t dimensions, const unsigned *bits_per_axis)
{
    size_t bits = curve_key_bits_mixed(dimensions, bits_per_axis);
    if (bits == 0) return 0;
    return (bits + 63u) / 64u;
}

CURVE_INLINE curve_status curve_mortonnd_encode_words_mixed(
    const uint32_t *coords, size_t dimensions,
    const unsigned *bits_per_axis,
    uint64_t *key_words, size_t key_capacity_words)
{
    curve_status st = curve__check_mixed_bits(dimensions, bits_per_axis);
    size_t need, i, k, W, K;
    unsigned maxb, p;
    if (st != CURVE_OK) return st;
    if (coords == NULL || key_words == NULL) return CURVE_EINVAL;
    K = curve_key_bits_mixed(dimensions, bits_per_axis);
    if (K == 0u) return CURVE_EOVERFLOW;
    need = (K + 63u) / 64u;
    if (key_capacity_words < need) return CURVE_ENOSPACE;
    W = need;
    for (i = 0; i < W; i++) key_words[i] = 0;
    maxb = 0;
    for (i = 0; i < dimensions; i++)
        if (bits_per_axis[i] > maxb) maxb = bits_per_axis[i];
    k = 0;
    p = maxb;
    while (p-- > 0u) {
        for (i = 0; i < dimensions; i++) {
            uint32_t mask;
            uint64_t bit;
            if (bits_per_axis[i] <= p) continue;
            mask = curve__lo_mask((int)bits_per_axis[i]);
            bit = ((uint64_t)(coords[i] & mask) >> p) & 1ull;
            key_words[k / 64u] |= bit << (63u - (unsigned)(k % 64u));
            k++;
        }
    }
    return CURVE_OK;
}

CURVE_INLINE curve_status curve_mortonnd_decode_words_mixed(
    const uint64_t *key_words, size_t key_count_words,
    size_t dimensions, const unsigned *bits_per_axis,
    uint32_t *coords, size_t coords_capacity)
{
    curve_status st = curve__check_mixed_bits(dimensions, bits_per_axis);
    size_t K, i, k;
    unsigned maxb, p;
    if (st != CURVE_OK) return st;
    if (key_words == NULL || coords == NULL) return CURVE_EINVAL;
    if (coords_capacity < dimensions) return CURVE_ENOSPACE;
    K = curve_key_bits_mixed(dimensions, bits_per_axis);
    if (K == 0u) return CURVE_EOVERFLOW;
    st = curve__key_check_canonical(key_words, key_count_words, K);
    if (st != CURVE_OK) return st;
    for (i = 0; i < dimensions; i++) coords[i] = 0;
    maxb = 0;
    for (i = 0; i < dimensions; i++)
        if (bits_per_axis[i] > maxb) maxb = bits_per_axis[i];
    k = 0;
    p = maxb;
    while (p-- > 0u) {
        for (i = 0; i < dimensions; i++) {
            uint32_t bit;
            if (bits_per_axis[i] <= p) continue;
            bit = (uint32_t)((key_words[k / 64u] >>
                (63u - (unsigned)(k % 64u))) & 1ull);
            coords[i] |= bit << p;
            k++;
        }
    }
    return CURVE_OK;
}

/* -------------------------------------------------------------------------- */
/* Linear quantizers (float / double → integer bins)                          */
/*                                                                            */
/* Maps each axis independently from [min, max] onto 0 .. 2^bits - 1.         */
/* bits is 1..32. No <math.h>. Compile without -ffast-math for cross-arch     */
/* determinism (IEEE comparisons and round-to-nearest ties).                  */
/*                                                                            */
/* CURVE_QUANT_FLOOR   — geometric / GIS bins. n_bins = 2^bits equal-width    */
/*                       intervals on [min, max]; max maps to the last bin.   */
/* CURVE_QUANT_NEAREST — embedding grid. Affine map of [min, max] onto        */
/*                       0 .. 2^bits-1, then floor(s + 0.5); ties toward +∞.  */
/*                                                                            */
/* Non-finite inputs: NaN and −Inf → 0; +Inf → last bin. Finite x below min   */
/* → 0; above max → last. Zero-width finite [min, min] → 0 for every x.       */
/* Bounds must be finite with max >= min; otherwise CURVE_EINVAL.             */
/*                                                                            */
/* n_lo / n_hi / n_nan, when non-NULL, are *set* to this call's totals (not   */
/* accumulated). Pass NULL to skip. Fused encode quantizes into caller        */
/* scratch then packs a wide Morton/Hilbert key; output is bit-identical to   */
/* quantize + encode_words. Hilbert still materializes every axis (Skilling). */
/* f32 APIs promote to f64 (exact for finite floats).                         */
/* -------------------------------------------------------------------------- */

typedef enum curve_quant_round {
    CURVE_QUANT_FLOOR   = 0,
    CURVE_QUANT_NEAREST = 1
} curve_quant_round;

CURVE_INLINE int curve__dbl_isnan(double x)
{
    return x != x;
}

CURVE_INLINE int curve__dbl_isfinite(double x)
{
    return x == x && (x * 2.0 != x || x == 0.0);
}

CURVE_INLINE curve_status curve__quant_check_mode(unsigned bits,
                                                  curve_quant_round rnd)
{
    if (bits < 1u || bits > CURVE_MAX_BITS_PER_AXIS) return CURVE_EINVAL;
    if (rnd != CURVE_QUANT_FLOOR && rnd != CURVE_QUANT_NEAREST)
        return CURVE_EINVAL;
    return CURVE_OK;
}

CURVE_INLINE curve_status curve__quant_check_bounds(double lo, double hi)
{
    if (!curve__dbl_isfinite(lo) || !curve__dbl_isfinite(hi))
        return CURVE_EINVAL;
    if (hi < lo) return CURVE_EINVAL;
    return CURVE_OK;
}

CURVE_INLINE void curve__quant_store_counts(size_t *n_lo, size_t *n_hi,
                                            size_t *n_nan,
                                            size_t c_lo, size_t c_hi,
                                            size_t c_nan)
{
    if (n_lo) *n_lo = c_lo;
    if (n_hi) *n_hi = c_hi;
    if (n_nan) *n_nan = c_nan;
}

/* Requires bits in 1..32, rnd in {FLOOR,NEAREST}, finite lo, hi with hi>=lo. */
CURVE_INLINE uint32_t curve__quantize1_f64_raw(
    double x, double lo, double hi, unsigned bits,
    curve_quant_round rnd,
    size_t *c_lo, size_t *c_hi, size_t *c_nan)
{
    uint32_t last = curve__lo_mask((int)bits);
    double t, s;
    uint32_t q;

    if (curve__dbl_isnan(x)) {
        (*c_nan)++;
        return 0u;
    }
    if (x < lo) {
        (*c_lo)++;
        return 0u;
    }
    if (x > hi) {
        (*c_hi)++;
        return (hi == lo) ? 0u : last;
    }
    if (hi == lo) return 0u;

    t = (x - lo) / (hi - lo);
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;

    if (rnd == CURVE_QUANT_FLOOR)
        s = t * (double)(1ull << bits);
    else
        s = t * (double)last + 0.5;

    /* 2^32 is exact in IEEE-754 binary64. Do not cast s >= 2^32 to uint32. */
    if (!(s < 4294967296.0)) return last;
    q = (uint32_t)s;
    return (q > last) ? last : q;
}

CURVE_INLINE double curve__dequantize1_f64_raw(
    uint32_t q, double lo, double hi, unsigned bits, curve_quant_round rnd)
{
    uint32_t last = curve__lo_mask((int)bits);
    if (q > last) q = last;
    if (hi == lo) return lo;
    if (rnd == CURVE_QUANT_FLOOR) {
        double n_bins = (double)(1ull << bits);
        return lo + ((double)q + 0.5) * (hi - lo) / n_bins;
    }
    return lo + ((double)q / (double)last) * (hi - lo);
}

CURVE_INLINE curve_status curve__quantize_axes_f64(
    const double *CURVE_RESTRICT values, size_t n,
    const double *CURVE_RESTRICT mins,
    const double *CURVE_RESTRICT maxs,
    unsigned bits, curve_quant_round rnd,
    uint32_t *CURVE_RESTRICT out,
    size_t *c_lo, size_t *c_hi, size_t *c_nan)
{
    size_t i;
    for (i = 0; i < n; i++) {
        if (curve__quant_check_bounds(mins[i], maxs[i]) != CURVE_OK)
            return CURVE_EINVAL;
    }
    *c_lo = 0;
    *c_hi = 0;
    *c_nan = 0;
    for (i = 0; i < n; i++) {
        out[i] = curve__quantize1_f64_raw(values[i], mins[i], maxs[i],
                                          bits, rnd, c_lo, c_hi, c_nan);
    }
    return CURVE_OK;
}

CURVE_INLINE curve_status curve__quantize_axes_f32(
    const float *CURVE_RESTRICT values, size_t n,
    const float *CURVE_RESTRICT mins,
    const float *CURVE_RESTRICT maxs,
    unsigned bits, curve_quant_round rnd,
    uint32_t *CURVE_RESTRICT out,
    size_t *c_lo, size_t *c_hi, size_t *c_nan)
{
    size_t i;
    for (i = 0; i < n; i++) {
        if (curve__quant_check_bounds((double)mins[i], (double)maxs[i]) !=
            CURVE_OK)
            return CURVE_EINVAL;
    }
    *c_lo = 0;
    *c_hi = 0;
    *c_nan = 0;
    for (i = 0; i < n; i++) {
        out[i] = curve__quantize1_f64_raw((double)values[i], (double)mins[i],
                                          (double)maxs[i], bits, rnd,
                                          c_lo, c_hi, c_nan);
    }
    return CURVE_OK;
}

CURVE_INLINE curve_status curve_quantize1_f64(
    double x, double lo, double hi, unsigned bits,
    curve_quant_round rnd, uint32_t *out)
{
    size_t c_lo = 0, c_hi = 0, c_nan = 0;
    curve_status st = curve__quant_check_mode(bits, rnd);
    if (st != CURVE_OK) return st;
    if (out == NULL) return CURVE_EINVAL;
    st = curve__quant_check_bounds(lo, hi);
    if (st != CURVE_OK) return st;
    *out = curve__quantize1_f64_raw(x, lo, hi, bits, rnd,
                                    &c_lo, &c_hi, &c_nan);
    return CURVE_OK;
}

CURVE_INLINE curve_status curve_quantize1_f32(
    float x, float lo, float hi, unsigned bits,
    curve_quant_round rnd, uint32_t *out)
{
    return curve_quantize1_f64((double)x, (double)lo, (double)hi,
                               bits, rnd, out);
}

CURVE_INLINE curve_status curve_quantize_f64(
    const double *CURVE_RESTRICT values, size_t n,
    const double *CURVE_RESTRICT mins,
    const double *CURVE_RESTRICT maxs,
    unsigned bits, curve_quant_round rnd,
    uint32_t *CURVE_RESTRICT out, size_t cap,
    size_t *n_lo, size_t *n_hi, size_t *n_nan)
{
    curve_status st;
    size_t c_lo = 0, c_hi = 0, c_nan = 0;
    st = curve__quant_check_mode(bits, rnd);
    if (st != CURVE_OK) return st;
    if (n == 0u) {
        curve__quant_store_counts(n_lo, n_hi, n_nan, 0, 0, 0);
        return CURVE_OK;
    }
    if (values == NULL || mins == NULL || maxs == NULL || out == NULL)
        return CURVE_EINVAL;
    if (cap < n) return CURVE_ENOSPACE;
    st = curve__quantize_axes_f64(values, n, mins, maxs, bits, rnd, out,
                                  &c_lo, &c_hi, &c_nan);
    if (st != CURVE_OK) return st;
    curve__quant_store_counts(n_lo, n_hi, n_nan, c_lo, c_hi, c_nan);
    return CURVE_OK;
}

CURVE_INLINE curve_status curve_quantize_f32(
    const float *CURVE_RESTRICT values, size_t n,
    const float *CURVE_RESTRICT mins,
    const float *CURVE_RESTRICT maxs,
    unsigned bits, curve_quant_round rnd,
    uint32_t *CURVE_RESTRICT out, size_t cap,
    size_t *n_lo, size_t *n_hi, size_t *n_nan)
{
    curve_status st;
    size_t c_lo = 0, c_hi = 0, c_nan = 0;
    st = curve__quant_check_mode(bits, rnd);
    if (st != CURVE_OK) return st;
    if (n == 0u) {
        curve__quant_store_counts(n_lo, n_hi, n_nan, 0, 0, 0);
        return CURVE_OK;
    }
    if (values == NULL || mins == NULL || maxs == NULL || out == NULL)
        return CURVE_EINVAL;
    if (cap < n) return CURVE_ENOSPACE;
    st = curve__quantize_axes_f32(values, n, mins, maxs, bits, rnd, out,
                                  &c_lo, &c_hi, &c_nan);
    if (st != CURVE_OK) return st;
    curve__quant_store_counts(n_lo, n_hi, n_nan, c_lo, c_hi, c_nan);
    return CURVE_OK;
}

CURVE_INLINE curve_status curve_dequantize1_f64(
    uint32_t q, double lo, double hi, unsigned bits,
    curve_quant_round rnd, double *out)
{
    curve_status st = curve__quant_check_mode(bits, rnd);
    uint32_t last;
    if (st != CURVE_OK) return st;
    if (out == NULL) return CURVE_EINVAL;
    st = curve__quant_check_bounds(lo, hi);
    if (st != CURVE_OK) return st;
    last = curve__lo_mask((int)bits);
    if (q > last) return CURVE_EINVAL;
    *out = curve__dequantize1_f64_raw(q, lo, hi, bits, rnd);
    return CURVE_OK;
}

CURVE_INLINE curve_status curve_dequantize1_f32(
    uint32_t q, float lo, float hi, unsigned bits,
    curve_quant_round rnd, float *out)
{
    double y;
    curve_status st;
    if (out == NULL) return CURVE_EINVAL;
    st = curve_dequantize1_f64(q, (double)lo, (double)hi, bits, rnd, &y);
    if (st != CURVE_OK) return st;
    *out = (float)y;
    return CURVE_OK;
}

CURVE_INLINE curve_status curve_dequantize_f64(
    const uint32_t *CURVE_RESTRICT qs, size_t n,
    const double *CURVE_RESTRICT mins,
    const double *CURVE_RESTRICT maxs,
    unsigned bits, curve_quant_round rnd,
    double *CURVE_RESTRICT out, size_t cap)
{
    curve_status st;
    uint32_t last;
    size_t i;
    st = curve__quant_check_mode(bits, rnd);
    if (st != CURVE_OK) return st;
    if (n == 0u) return CURVE_OK;
    if (qs == NULL || mins == NULL || maxs == NULL || out == NULL)
        return CURVE_EINVAL;
    if (cap < n) return CURVE_ENOSPACE;
    last = curve__lo_mask((int)bits);
    for (i = 0; i < n; i++) {
        if (curve__quant_check_bounds(mins[i], maxs[i]) != CURVE_OK)
            return CURVE_EINVAL;
        if (qs[i] > last) return CURVE_EINVAL;
    }
    for (i = 0; i < n; i++)
        out[i] = curve__dequantize1_f64_raw(qs[i], mins[i], maxs[i],
                                            bits, rnd);
    return CURVE_OK;
}

CURVE_INLINE curve_status curve_dequantize_f32(
    const uint32_t *CURVE_RESTRICT qs, size_t n,
    const float *CURVE_RESTRICT mins,
    const float *CURVE_RESTRICT maxs,
    unsigned bits, curve_quant_round rnd,
    float *CURVE_RESTRICT out, size_t cap)
{
    curve_status st;
    uint32_t last;
    size_t i;
    st = curve__quant_check_mode(bits, rnd);
    if (st != CURVE_OK) return st;
    if (n == 0u) return CURVE_OK;
    if (qs == NULL || mins == NULL || maxs == NULL || out == NULL)
        return CURVE_EINVAL;
    if (cap < n) return CURVE_ENOSPACE;
    last = curve__lo_mask((int)bits);
    for (i = 0; i < n; i++) {
        if (curve__quant_check_bounds((double)mins[i], (double)maxs[i]) !=
            CURVE_OK)
            return CURVE_EINVAL;
        if (qs[i] > last) return CURVE_EINVAL;
    }
    for (i = 0; i < n; i++)
        out[i] = (float)curve__dequantize1_f64_raw(
            qs[i], (double)mins[i], (double)maxs[i], bits, rnd);
    return CURVE_OK;
}

CURVE_INLINE curve_status curve_mortonnd_encode_f64(
    const double *CURVE_RESTRICT values, size_t dimensions,
    const double *CURVE_RESTRICT mins,
    const double *CURVE_RESTRICT maxs,
    unsigned bits_per_axis, curve_quant_round rnd,
    uint64_t *CURVE_RESTRICT key_words, size_t key_capacity_words,
    uint32_t *CURVE_RESTRICT scratch_axes, size_t scratch_capacity_axes,
    size_t *n_lo, size_t *n_hi, size_t *n_nan)
{
    curve_status st;
    size_t c_lo = 0, c_hi = 0, c_nan = 0, need;
    st = curve__check_dims_bits(dimensions, bits_per_axis);
    if (st != CURVE_OK) return st;
    st = curve__quant_check_mode(bits_per_axis, rnd);
    if (st != CURVE_OK) return st;
    if (values == NULL || mins == NULL || maxs == NULL || key_words == NULL ||
        scratch_axes == NULL)
        return CURVE_EINVAL;
    need = curve_key_words_required(dimensions, bits_per_axis);
    if (need == 0u) return CURVE_EOVERFLOW;
    if (key_capacity_words < need) return CURVE_ENOSPACE;
    if (scratch_capacity_axes < dimensions) return CURVE_ENOSPACE;
    st = curve__quantize_axes_f64(values, dimensions, mins, maxs,
                                  bits_per_axis, rnd, scratch_axes,
                                  &c_lo, &c_hi, &c_nan);
    if (st != CURVE_OK) return st;
    st = curve_mortonnd_encode_words(scratch_axes, dimensions, bits_per_axis,
                                     key_words, key_capacity_words);
    if (st != CURVE_OK) return st;
    curve__quant_store_counts(n_lo, n_hi, n_nan, c_lo, c_hi, c_nan);
    return CURVE_OK;
}

CURVE_INLINE curve_status curve_hilbertnd_encode_f64(
    const double *CURVE_RESTRICT values, size_t dimensions,
    const double *CURVE_RESTRICT mins,
    const double *CURVE_RESTRICT maxs,
    unsigned bits_per_axis, curve_quant_round rnd,
    uint64_t *CURVE_RESTRICT key_words, size_t key_capacity_words,
    uint32_t *CURVE_RESTRICT scratch_axes, size_t scratch_capacity_axes,
    size_t *n_lo, size_t *n_hi, size_t *n_nan)
{
    curve_status st;
    size_t c_lo = 0, c_hi = 0, c_nan = 0, need;
    st = curve__check_dims_bits(dimensions, bits_per_axis);
    if (st != CURVE_OK) return st;
    st = curve__quant_check_mode(bits_per_axis, rnd);
    if (st != CURVE_OK) return st;
    if (values == NULL || mins == NULL || maxs == NULL || key_words == NULL ||
        scratch_axes == NULL)
        return CURVE_EINVAL;
    need = curve_key_words_required(dimensions, bits_per_axis);
    if (need == 0u) return CURVE_EOVERFLOW;
    if (key_capacity_words < need) return CURVE_ENOSPACE;
    if (scratch_capacity_axes < dimensions) return CURVE_ENOSPACE;
    st = curve__quantize_axes_f64(values, dimensions, mins, maxs,
                                  bits_per_axis, rnd, scratch_axes,
                                  &c_lo, &c_hi, &c_nan);
    if (st != CURVE_OK) return st;
    st = curve_hilbertnd_encode_words(scratch_axes, dimensions, bits_per_axis,
                                      key_words, key_capacity_words,
                                      scratch_axes, scratch_capacity_axes);
    if (st != CURVE_OK) return st;
    curve__quant_store_counts(n_lo, n_hi, n_nan, c_lo, c_hi, c_nan);
    return CURVE_OK;
}

CURVE_INLINE curve_status curve_mortonnd_encode_f32(
    const float *CURVE_RESTRICT values, size_t dimensions,
    const float *CURVE_RESTRICT mins,
    const float *CURVE_RESTRICT maxs,
    unsigned bits_per_axis, curve_quant_round rnd,
    uint64_t *CURVE_RESTRICT key_words, size_t key_capacity_words,
    uint32_t *CURVE_RESTRICT scratch_axes, size_t scratch_capacity_axes,
    size_t *n_lo, size_t *n_hi, size_t *n_nan)
{
    curve_status st;
    size_t c_lo = 0, c_hi = 0, c_nan = 0, need;
    st = curve__check_dims_bits(dimensions, bits_per_axis);
    if (st != CURVE_OK) return st;
    st = curve__quant_check_mode(bits_per_axis, rnd);
    if (st != CURVE_OK) return st;
    if (values == NULL || mins == NULL || maxs == NULL || key_words == NULL ||
        scratch_axes == NULL)
        return CURVE_EINVAL;
    need = curve_key_words_required(dimensions, bits_per_axis);
    if (need == 0u) return CURVE_EOVERFLOW;
    if (key_capacity_words < need) return CURVE_ENOSPACE;
    if (scratch_capacity_axes < dimensions) return CURVE_ENOSPACE;
    st = curve__quantize_axes_f32(values, dimensions, mins, maxs,
                                  bits_per_axis, rnd, scratch_axes,
                                  &c_lo, &c_hi, &c_nan);
    if (st != CURVE_OK) return st;
    st = curve_mortonnd_encode_words(scratch_axes, dimensions, bits_per_axis,
                                     key_words, key_capacity_words);
    if (st != CURVE_OK) return st;
    curve__quant_store_counts(n_lo, n_hi, n_nan, c_lo, c_hi, c_nan);
    return CURVE_OK;
}

CURVE_INLINE curve_status curve_hilbertnd_encode_f32(
    const float *CURVE_RESTRICT values, size_t dimensions,
    const float *CURVE_RESTRICT mins,
    const float *CURVE_RESTRICT maxs,
    unsigned bits_per_axis, curve_quant_round rnd,
    uint64_t *CURVE_RESTRICT key_words, size_t key_capacity_words,
    uint32_t *CURVE_RESTRICT scratch_axes, size_t scratch_capacity_axes,
    size_t *n_lo, size_t *n_hi, size_t *n_nan)
{
    curve_status st;
    size_t c_lo = 0, c_hi = 0, c_nan = 0, need;
    st = curve__check_dims_bits(dimensions, bits_per_axis);
    if (st != CURVE_OK) return st;
    st = curve__quant_check_mode(bits_per_axis, rnd);
    if (st != CURVE_OK) return st;
    if (values == NULL || mins == NULL || maxs == NULL || key_words == NULL ||
        scratch_axes == NULL)
        return CURVE_EINVAL;
    need = curve_key_words_required(dimensions, bits_per_axis);
    if (need == 0u) return CURVE_EOVERFLOW;
    if (key_capacity_words < need) return CURVE_ENOSPACE;
    if (scratch_capacity_axes < dimensions) return CURVE_ENOSPACE;
    st = curve__quantize_axes_f32(values, dimensions, mins, maxs,
                                  bits_per_axis, rnd, scratch_axes,
                                  &c_lo, &c_hi, &c_nan);
    if (st != CURVE_OK) return st;
    st = curve_hilbertnd_encode_words(scratch_axes, dimensions, bits_per_axis,
                                      key_words, key_capacity_words,
                                      scratch_axes, scratch_capacity_axes);
    if (st != CURVE_OK) return st;
    curve__quant_store_counts(n_lo, n_hi, n_nan, c_lo, c_hi, c_nan);
    return CURVE_OK;
}

/* AoS batch: values[i * value_stride + d], keys[i * key_stride_words + w].   */
/* Shared per-axis mins/maxs. Scratch holds one vector (dimensions uint32).   */
CURVE_INLINE curve_status curve_mortonnd_encode_f64_aos(
    const double *CURVE_RESTRICT values, size_t n_vecs, size_t value_stride,
    size_t dimensions,
    const double *CURVE_RESTRICT mins,
    const double *CURVE_RESTRICT maxs,
    unsigned bits_per_axis, curve_quant_round rnd,
    uint64_t *CURVE_RESTRICT keys, size_t key_stride_words,
    uint32_t *CURVE_RESTRICT scratch_axes, size_t scratch_capacity_axes,
    size_t *n_lo, size_t *n_hi, size_t *n_nan)
{
    curve_status st;
    size_t i, need, t_lo = 0, t_hi = 0, t_nan = 0;
    st = curve__check_dims_bits(dimensions, bits_per_axis);
    if (st != CURVE_OK) return st;
    st = curve__quant_check_mode(bits_per_axis, rnd);
    if (st != CURVE_OK) return st;
    if (n_vecs == 0u) {
        curve__quant_store_counts(n_lo, n_hi, n_nan, 0, 0, 0);
        return CURVE_OK;
    }
    if (values == NULL || mins == NULL || maxs == NULL || keys == NULL ||
        scratch_axes == NULL)
        return CURVE_EINVAL;
    if (value_stride < dimensions) return CURVE_EINVAL;
    need = curve_key_words_required(dimensions, bits_per_axis);
    if (need == 0u) return CURVE_EOVERFLOW;
    if (key_stride_words < need) return CURVE_EINVAL;
    if (scratch_capacity_axes < dimensions) return CURVE_ENOSPACE;
    if (value_stride != 0u && n_vecs > ((size_t)-1) / value_stride)
        return CURVE_EOVERFLOW;
    if (key_stride_words != 0u && n_vecs > ((size_t)-1) / key_stride_words)
        return CURVE_EOVERFLOW;
    for (i = 0; i < n_vecs; i++) {
        size_t c_lo = 0, c_hi = 0, c_nan = 0;
        st = curve__quantize_axes_f64(values + i * value_stride, dimensions,
                                      mins, maxs, bits_per_axis, rnd,
                                      scratch_axes, &c_lo, &c_hi, &c_nan);
        if (st != CURVE_OK) return st;
        st = curve_mortonnd_encode_words(scratch_axes, dimensions,
                                         bits_per_axis,
                                         keys + i * key_stride_words,
                                         key_stride_words);
        if (st != CURVE_OK) return st;
        t_lo += c_lo;
        t_hi += c_hi;
        t_nan += c_nan;
    }
    curve__quant_store_counts(n_lo, n_hi, n_nan, t_lo, t_hi, t_nan);
    return CURVE_OK;
}

CURVE_INLINE curve_status curve_hilbertnd_encode_f64_aos(
    const double *CURVE_RESTRICT values, size_t n_vecs, size_t value_stride,
    size_t dimensions,
    const double *CURVE_RESTRICT mins,
    const double *CURVE_RESTRICT maxs,
    unsigned bits_per_axis, curve_quant_round rnd,
    uint64_t *CURVE_RESTRICT keys, size_t key_stride_words,
    uint32_t *CURVE_RESTRICT scratch_axes, size_t scratch_capacity_axes,
    size_t *n_lo, size_t *n_hi, size_t *n_nan)
{
    curve_status st;
    size_t i, need, t_lo = 0, t_hi = 0, t_nan = 0;
    st = curve__check_dims_bits(dimensions, bits_per_axis);
    if (st != CURVE_OK) return st;
    st = curve__quant_check_mode(bits_per_axis, rnd);
    if (st != CURVE_OK) return st;
    if (n_vecs == 0u) {
        curve__quant_store_counts(n_lo, n_hi, n_nan, 0, 0, 0);
        return CURVE_OK;
    }
    if (values == NULL || mins == NULL || maxs == NULL || keys == NULL ||
        scratch_axes == NULL)
        return CURVE_EINVAL;
    if (value_stride < dimensions) return CURVE_EINVAL;
    need = curve_key_words_required(dimensions, bits_per_axis);
    if (need == 0u) return CURVE_EOVERFLOW;
    if (key_stride_words < need) return CURVE_EINVAL;
    if (scratch_capacity_axes < dimensions) return CURVE_ENOSPACE;
    if (value_stride != 0u && n_vecs > ((size_t)-1) / value_stride)
        return CURVE_EOVERFLOW;
    if (key_stride_words != 0u && n_vecs > ((size_t)-1) / key_stride_words)
        return CURVE_EOVERFLOW;
    for (i = 0; i < n_vecs; i++) {
        size_t c_lo = 0, c_hi = 0, c_nan = 0;
        st = curve__quantize_axes_f64(values + i * value_stride, dimensions,
                                      mins, maxs, bits_per_axis, rnd,
                                      scratch_axes, &c_lo, &c_hi, &c_nan);
        if (st != CURVE_OK) return st;
        st = curve_hilbertnd_encode_words(scratch_axes, dimensions,
                                          bits_per_axis,
                                          keys + i * key_stride_words,
                                          key_stride_words,
                                          scratch_axes, scratch_capacity_axes);
        if (st != CURVE_OK) return st;
        t_lo += c_lo;
        t_hi += c_hi;
        t_nan += c_nan;
    }
    curve__quant_store_counts(n_lo, n_hi, n_nan, t_lo, t_hi, t_nan);
    return CURVE_OK;
}

CURVE_INLINE curve_status curve_mortonnd_encode_f32_aos(
    const float *CURVE_RESTRICT values, size_t n_vecs, size_t value_stride,
    size_t dimensions,
    const float *CURVE_RESTRICT mins,
    const float *CURVE_RESTRICT maxs,
    unsigned bits_per_axis, curve_quant_round rnd,
    uint64_t *CURVE_RESTRICT keys, size_t key_stride_words,
    uint32_t *CURVE_RESTRICT scratch_axes, size_t scratch_capacity_axes,
    size_t *n_lo, size_t *n_hi, size_t *n_nan)
{
    curve_status st;
    size_t i, need, t_lo = 0, t_hi = 0, t_nan = 0;
    st = curve__check_dims_bits(dimensions, bits_per_axis);
    if (st != CURVE_OK) return st;
    st = curve__quant_check_mode(bits_per_axis, rnd);
    if (st != CURVE_OK) return st;
    if (n_vecs == 0u) {
        curve__quant_store_counts(n_lo, n_hi, n_nan, 0, 0, 0);
        return CURVE_OK;
    }
    if (values == NULL || mins == NULL || maxs == NULL || keys == NULL ||
        scratch_axes == NULL)
        return CURVE_EINVAL;
    if (value_stride < dimensions) return CURVE_EINVAL;
    need = curve_key_words_required(dimensions, bits_per_axis);
    if (need == 0u) return CURVE_EOVERFLOW;
    if (key_stride_words < need) return CURVE_EINVAL;
    if (scratch_capacity_axes < dimensions) return CURVE_ENOSPACE;
    if (value_stride != 0u && n_vecs > ((size_t)-1) / value_stride)
        return CURVE_EOVERFLOW;
    if (key_stride_words != 0u && n_vecs > ((size_t)-1) / key_stride_words)
        return CURVE_EOVERFLOW;
    for (i = 0; i < n_vecs; i++) {
        size_t c_lo = 0, c_hi = 0, c_nan = 0;
        st = curve__quantize_axes_f32(values + i * value_stride, dimensions,
                                      mins, maxs, bits_per_axis, rnd,
                                      scratch_axes, &c_lo, &c_hi, &c_nan);
        if (st != CURVE_OK) return st;
        st = curve_mortonnd_encode_words(scratch_axes, dimensions,
                                         bits_per_axis,
                                         keys + i * key_stride_words,
                                         key_stride_words);
        if (st != CURVE_OK) return st;
        t_lo += c_lo;
        t_hi += c_hi;
        t_nan += c_nan;
    }
    curve__quant_store_counts(n_lo, n_hi, n_nan, t_lo, t_hi, t_nan);
    return CURVE_OK;
}

CURVE_INLINE curve_status curve_hilbertnd_encode_f32_aos(
    const float *CURVE_RESTRICT values, size_t n_vecs, size_t value_stride,
    size_t dimensions,
    const float *CURVE_RESTRICT mins,
    const float *CURVE_RESTRICT maxs,
    unsigned bits_per_axis, curve_quant_round rnd,
    uint64_t *CURVE_RESTRICT keys, size_t key_stride_words,
    uint32_t *CURVE_RESTRICT scratch_axes, size_t scratch_capacity_axes,
    size_t *n_lo, size_t *n_hi, size_t *n_nan)
{
    curve_status st;
    size_t i, need, t_lo = 0, t_hi = 0, t_nan = 0;
    st = curve__check_dims_bits(dimensions, bits_per_axis);
    if (st != CURVE_OK) return st;
    st = curve__quant_check_mode(bits_per_axis, rnd);
    if (st != CURVE_OK) return st;
    if (n_vecs == 0u) {
        curve__quant_store_counts(n_lo, n_hi, n_nan, 0, 0, 0);
        return CURVE_OK;
    }
    if (values == NULL || mins == NULL || maxs == NULL || keys == NULL ||
        scratch_axes == NULL)
        return CURVE_EINVAL;
    if (value_stride < dimensions) return CURVE_EINVAL;
    need = curve_key_words_required(dimensions, bits_per_axis);
    if (need == 0u) return CURVE_EOVERFLOW;
    if (key_stride_words < need) return CURVE_EINVAL;
    if (scratch_capacity_axes < dimensions) return CURVE_ENOSPACE;
    if (value_stride != 0u && n_vecs > ((size_t)-1) / value_stride)
        return CURVE_EOVERFLOW;
    if (key_stride_words != 0u && n_vecs > ((size_t)-1) / key_stride_words)
        return CURVE_EOVERFLOW;
    for (i = 0; i < n_vecs; i++) {
        size_t c_lo = 0, c_hi = 0, c_nan = 0;
        st = curve__quantize_axes_f32(values + i * value_stride, dimensions,
                                      mins, maxs, bits_per_axis, rnd,
                                      scratch_axes, &c_lo, &c_hi, &c_nan);
        if (st != CURVE_OK) return st;
        st = curve_hilbertnd_encode_words(scratch_axes, dimensions,
                                          bits_per_axis,
                                          keys + i * key_stride_words,
                                          key_stride_words,
                                          scratch_axes, scratch_capacity_axes);
        if (st != CURVE_OK) return st;
        t_lo += c_lo;
        t_hi += c_hi;
        t_nan += c_nan;
    }
    curve__quant_store_counts(n_lo, n_hi, n_nan, t_lo, t_hi, t_nan);
    return CURVE_OK;
}

/* -------------------------------------------------------------------------- */
/* Spatial cells — complete-level prefixes of a wide key                      */
/*                                                                            */
/* A cell is a prefix whose length is a multiple of `dimensions` (one Hilbert */
/* / Morton level = n_dims bits). Partial-level prefixes are rejected with    */
/* CURVE_EINVAL. Prefix words are caller-owned.                               */
/*                                                                            */
/* Child indices are in key order: child 0 is the first sub-cell along the    */
/* curve. `curve_cell_child` / sibling / child_count require n_dims <= 64     */
/* (the child index is a uint64_t). Bounds APIs work at any dimension count.  */
/* Morton bounds are read from interleaved prefix bits. Hilbert bounds decode */
/* one point in the cell and snap to the dyadic hypercube.                    */
/* -------------------------------------------------------------------------- */

typedef struct curve_cell {
    size_t   dimensions;
    unsigned bits_per_axis;
    size_t   prefix_bits; /* multiple of dimensions, 0 = root */
} curve_cell;

CURVE_INLINE curve_status curve_cell_init(curve_cell *cell, size_t dimensions,
                                          unsigned bits_per_axis,
                                          size_t prefix_bits)
{
    curve_status st;
    size_t max_bits;
    if (cell == NULL) return CURVE_EINVAL;
    st = curve__check_dims_bits(dimensions, bits_per_axis);
    if (st != CURVE_OK) return st;
    if (prefix_bits % dimensions != 0u) return CURVE_EINVAL;
    max_bits = dimensions * (size_t)bits_per_axis;
    if (prefix_bits > max_bits) return CURVE_EINVAL;
    cell->dimensions = dimensions;
    cell->bits_per_axis = bits_per_axis;
    cell->prefix_bits = prefix_bits;
    return CURVE_OK;
}

CURVE_INLINE size_t curve_cell_depth(const curve_cell *cell)
{
    if (cell == NULL || cell->dimensions == 0u) return 0;
    return cell->prefix_bits / cell->dimensions;
}

CURVE_INLINE size_t curve_cell_max_depth(const curve_cell *cell)
{
    if (cell == NULL) return 0;
    return (size_t)cell->bits_per_axis;
}

CURVE_INLINE int curve_cell_is_leaf(const curve_cell *cell)
{
    if (cell == NULL) return 0;
    return curve_cell_depth(cell) == (size_t)cell->bits_per_axis;
}

CURVE_INLINE curve_status curve_cell_child_count(const curve_cell *cell,
                                                 size_t *out_count)
{
    size_t n;
    if (cell == NULL || out_count == NULL) return CURVE_EINVAL;
    if (curve_cell_is_leaf(cell)) {
        *out_count = 0;
        return CURVE_OK;
    }
    n = cell->dimensions;
    if (n > 64u) return CURVE_EUNSUPPORTED;
    if (n >= sizeof(size_t) * 8u) return CURVE_EOVERFLOW;
    *out_count = (size_t)1 << n;
    return CURVE_OK;
}

CURVE_INLINE curve_status curve__cell_require(const curve_cell *cell,
                                              const uint64_t *prefix)
{
    size_t W, key_bits;
    if (cell == NULL || prefix == NULL) return CURVE_EINVAL;
    if (cell->dimensions < 1u) return CURVE_EINVAL;
    if (cell->prefix_bits % cell->dimensions != 0u) return CURVE_EINVAL;
    key_bits = cell->dimensions * (size_t)cell->bits_per_axis;
    W = curve__key_n_words(key_bits);
    return curve__key_check_canonical(prefix, W, key_bits);
}

CURVE_INLINE curve_status curve_cell_parent(const curve_cell *cell,
                                            const uint64_t *prefix,
                                            curve_cell *parent,
                                            uint64_t *parent_prefix,
                                            size_t parent_words)
{
    curve_status st;
    size_t key_bits;
    if (cell == NULL || parent == NULL) return CURVE_EINVAL;
    if (cell->prefix_bits < cell->dimensions) return CURVE_ERANGE;
    st = curve__cell_require(cell, prefix);
    if (st != CURVE_OK) return st;
    key_bits = cell->dimensions * (size_t)cell->bits_per_axis;
    st = curve_key_prefix_mask(prefix, key_bits,
                               cell->prefix_bits - cell->dimensions,
                               parent_prefix, parent_words);
    if (st != CURVE_OK) return st;
    return curve_cell_init(parent, cell->dimensions, cell->bits_per_axis,
                           cell->prefix_bits - cell->dimensions);
}

CURVE_INLINE curve_status curve_cell_child(const curve_cell *cell,
                                           const uint64_t *prefix,
                                           uint64_t child_index,
                                           curve_cell *child,
                                           uint64_t *child_prefix,
                                           size_t child_words)
{
    curve_status st;
    size_t key_bits, n, i, new_bits;
    uint64_t nchild;
    if (cell == NULL || child == NULL) return CURVE_EINVAL;
    if (curve_cell_is_leaf(cell)) return CURVE_ERANGE;
    n = cell->dimensions;
    if (n > 64u) return CURVE_EUNSUPPORTED;
    nchild = (n == 64u) ? ~0ull : ((1ull << n) - 1ull);
    if (n < 64u && child_index > nchild) return CURVE_EINVAL;
    st = curve__cell_require(cell, prefix);
    if (st != CURVE_OK) return st;
    key_bits = n * (size_t)cell->bits_per_axis;
    new_bits = cell->prefix_bits + n;
    st = curve_key_prefix_mask(prefix, key_bits, cell->prefix_bits,
                               child_prefix, child_words);
    if (st != CURVE_OK) return st;
    for (i = 0; i < n; i++) {
        unsigned bit = (unsigned)((child_index >> (n - 1u - i)) & 1ull);
        curve__key_set_bit(child_prefix, cell->prefix_bits + i, bit);
    }
    return curve_cell_init(child, cell->dimensions, cell->bits_per_axis,
                           new_bits);
}

CURVE_INLINE curve_status curve_cell_child_index(const curve_cell *cell,
                                                 const uint64_t *prefix,
                                                 uint64_t *out_index)
{
    size_t n, i, start;
    uint64_t idx = 0;
    curve_status st;
    if (out_index == NULL) return CURVE_EINVAL;
    st = curve__cell_require(cell, prefix);
    if (st != CURVE_OK) return st;
    if (cell->prefix_bits < cell->dimensions) return CURVE_ERANGE;
    n = cell->dimensions;
    if (n > 64u) return CURVE_EUNSUPPORTED;
    start = cell->prefix_bits - n;
    for (i = 0; i < n; i++)
        idx = (idx << 1) | (uint64_t)curve__key_get_bit(prefix, start + i);
    *out_index = idx;
    return CURVE_OK;
}

CURVE_INLINE curve_status curve_cell_sibling(const curve_cell *cell,
                                             const uint64_t *prefix,
                                             int delta,
                                             curve_cell *sib,
                                             uint64_t *sib_prefix,
                                             size_t sib_words)
{
    curve_status st;
    uint64_t idx, nchild;
    size_t n, i, key_bits, W;
    curve_cell parent;
    if (delta == 0) {
        st = curve__cell_require(cell, prefix);
        if (st != CURVE_OK) return st;
        if (sib == NULL || sib_prefix == NULL) return CURVE_EINVAL;
        key_bits = cell->dimensions * (size_t)cell->bits_per_axis;
        W = curve__key_n_words(key_bits);
        if (sib_words < W) return CURVE_ENOSPACE;
        for (i = 0; i < W; i++) sib_prefix[i] = prefix[i];
        *sib = *cell;
        return CURVE_OK;
    }
    st = curve_cell_child_index(cell, prefix, &idx);
    if (st != CURVE_OK) return st;
    if (sib == NULL || sib_prefix == NULL) return CURVE_EINVAL;
    n = cell->dimensions;
    if (n >= 64u) return CURVE_EUNSUPPORTED;
    nchild = 1ull << n;
    if (delta > 0) {
        uint64_t d = (uint64_t)delta;
        if (d > nchild - 1ull - idx) return CURVE_ERANGE;
        idx += d;
    } else {
        uint64_t d = 0ull - (uint64_t)delta; /* |delta|, defined at INT_MIN */
        if (d > idx) return CURVE_ERANGE;
        idx -= d;
    }
    key_bits = n * (size_t)cell->bits_per_axis;
    W = curve__key_n_words(key_bits);
    if (sib_words < W) return CURVE_ENOSPACE;
    st = curve_cell_parent(cell, prefix, &parent, sib_prefix, sib_words);
    if (st != CURVE_OK) return st;
    return curve_cell_child(&parent, sib_prefix, idx, sib, sib_prefix, sib_words);
}

CURVE_INLINE void curve__cell_snap_bounds(const uint32_t *coords, size_t n,
                                          unsigned bits, size_t level,
                                          uint32_t *mins, uint32_t *maxs)
{
    uint32_t mask;
    size_t i;
    unsigned free_bits;
    free_bits = (level >= (size_t)bits) ? 0u : (unsigned)((size_t)bits - level);
    mask = curve__lo_mask((int)free_bits);
    for (i = 0; i < n; i++) {
        mins[i] = coords[i] & ~mask;
        maxs[i] = mins[i] | mask;
    }
}

CURVE_INLINE curve_status curve_morton_cell_bounds(const curve_cell *cell,
                                                   const uint64_t *prefix,
                                                   uint32_t *mins,
                                                   uint32_t *maxs,
                                                   size_t n_axes)
{
    curve_status st;
    size_t n, i, l, level;
    unsigned bits;
    uint32_t mask;
    if (mins == NULL || maxs == NULL) return CURVE_EINVAL;
    st = curve__cell_require(cell, prefix);
    if (st != CURVE_OK) return st;
    n = cell->dimensions;
    if (n_axes < n) return CURVE_ENOSPACE;
    bits = cell->bits_per_axis;
    level = curve_cell_depth(cell);
    for (i = 0; i < n; i++) mins[i] = maxs[i] = 0;
    for (l = 0; l < level; l++) {
        for (i = 0; i < n; i++) {
            unsigned bit = curve__key_get_bit(prefix, l * n + i);
            unsigned pos = (unsigned)(bits - 1u - (unsigned)l);
            uint32_t v = (uint32_t)bit << pos;
            mins[i] |= v;
            maxs[i] |= v;
        }
    }
    mask = curve__lo_mask((int)((level >= (size_t)bits)
                                    ? 0 : (int)((size_t)bits - level)));
    for (i = 0; i < n; i++) maxs[i] |= mask;
    return CURVE_OK;
}

CURVE_INLINE curve_status curve_hilbert_cell_bounds(const curve_cell *cell,
                                                    const uint64_t *prefix,
                                                    uint32_t *mins,
                                                    uint32_t *maxs,
                                                    size_t n_axes,
                                                    uint32_t *scratch,
                                                    size_t scratch_axes)
{
    curve_status st;
    size_t n, W, key_bits;
    if (mins == NULL || maxs == NULL) return CURVE_EINVAL;
    st = curve__cell_require(cell, prefix);
    if (st != CURVE_OK) return st;
    n = cell->dimensions;
    if (n_axes < n) return CURVE_ENOSPACE;
    if (scratch == NULL || scratch_axes < n) return CURVE_ENOSPACE;
    key_bits = n * (size_t)cell->bits_per_axis;
    W = curve__key_n_words(key_bits);
    /* Any point in the cell decodes to a coordinate inside the same
     * axis-aligned dyadic hypercube; snap drops bits below `depth`. */
    st = curve_hilbertnd_decode_words(prefix, W, n, cell->bits_per_axis,
                                      scratch, scratch_axes,
                                      scratch, scratch_axes);
    if (st != CURVE_OK) return st;
    curve__cell_snap_bounds(scratch, n, cell->bits_per_axis,
                            curve_cell_depth(cell), mins, maxs);
    return CURVE_OK;
}

/* -------------------------------------------------------------------------- */
/* AABB → key intervals (Morton / Hilbert)                                    */
/*                                                                            */
/* Inclusive integer box [box_min, box_max] on the 2^bits grid. Traversal     */
/* visits cells in key order, emits fully covered cells as prefix intervals   */
/* [lower, upper) (same convention as curve_key_prefix_interval), and merges  */
/* adjacent spans. n_dims is 1..CURVE_RANGE_MAX_DIMS (default 8).             */
/*                                                                            */
/* interval_capacity == 0 counts without writing (out_* may be NULL).         */
/* CURVE_RANGE_EXACT: if more intervals than capacity, returns ENOSPACE and   */
/* writes the first `capacity` spans; *out_count is the exact total.          */
/* CURVE_RANGE_OVERAPPROX: at most `capacity` spans that cover every in-box   */
/* point (false positives allowed). Never drops coverage.                     */
/* Scratch: curve_box_intervals_scratch_words / _scratch_axes.                */
/* -------------------------------------------------------------------------- */

typedef enum curve_range_mode {
    CURVE_RANGE_EXACT      = 0,
    CURVE_RANGE_OVERAPPROX = 1
} curve_range_mode;

CURVE_INLINE size_t curve_box_intervals_scratch_words(size_t dimensions,
                                                      unsigned bits_per_axis)
{
    size_t W, d, a, b;
    if (curve__check_dims_bits(dimensions, bits_per_axis) != CURVE_OK)
        return 0;
    W = curve_key_words_required(dimensions, bits_per_axis);
    d = (size_t)bits_per_axis + 1u;
    if (W > ((size_t)-1) / 4u) return 0;
    a = 4u * W;
    if ((W + 2u) > ((size_t)-1) / d) return 0;
    b = d * (W + 2u);
    if (a > ((size_t)-1) - b) return 0;
    return a + b;
}

CURVE_INLINE size_t curve_box_intervals_scratch_axes(size_t dimensions)
{
    if (dimensions < 1u || dimensions > (size_t)CURVE_MAX_DIMS)
        return 0;
    if (dimensions > ((size_t)-1) / 5u) return 0;
    return 5u * dimensions;
}

#define CURVE__AABB_DISJOINT  0
#define CURVE__AABB_CONTAINED 1
#define CURVE__AABB_PARTIAL   2

CURVE_INLINE int curve__aabb_relate(const uint32_t *cmins,
                                    const uint32_t *cmaxs,
                                    const uint32_t *bmins,
                                    const uint32_t *bmaxs, size_t n)
{
    size_t i;
    int contained = 1;
    for (i = 0; i < n; i++) {
        if (cmaxs[i] < bmins[i] || cmins[i] > bmaxs[i])
            return CURVE__AABB_DISJOINT;
        if (cmins[i] < bmins[i] || cmaxs[i] > bmaxs[i])
            contained = 0;
    }
    return contained ? CURVE__AABB_CONTAINED : CURVE__AABB_PARTIAL;
}

CURVE_INLINE curve_status curve__range_emit(
    const uint64_t *prefix, size_t key_bits, size_t prefix_bits, size_t W,
    uint64_t *last_lo, uint64_t *last_hi, int *last_unb, int *has_last,
    uint64_t *out_lo, uint64_t *out_hi, int *out_unb,
    size_t cap, size_t *written, size_t *count, curve_range_mode mode,
    uint64_t *ws_lo, uint64_t *ws_hi)
{
    int unb = 0, merge = 0;
    curve_status st;
    size_t i;
    st = curve_key_prefix_interval(prefix, key_bits, prefix_bits,
                                   ws_lo, ws_hi, W, &unb);
    if (st != CURVE_OK) return st;
    if (*has_last && *last_unb) return CURVE_OK;
    if (*has_last && curve_key_compare(ws_lo, last_hi, W) == 0) merge = 1;
    if (merge) {
        for (i = 0; i < W; i++) last_hi[i] = ws_hi[i];
        *last_unb = unb;
        if (cap > 0u && *written > 0u) {
            size_t slot = *written - 1u;
            for (i = 0; i < W; i++) out_hi[slot * W + i] = last_hi[i];
            out_unb[slot] = unb;
        }
        return CURVE_OK;
    }
    if (mode == CURVE_RANGE_OVERAPPROX && cap > 0u && *written >= cap &&
        *has_last) {
        for (i = 0; i < W; i++) last_hi[i] = ws_hi[i];
        *last_unb = unb;
        for (i = 0; i < W; i++)
            out_hi[(*written - 1u) * W + i] = ws_hi[i];
        out_unb[*written - 1u] = unb;
        return CURVE_OK;
    }
    if (cap > 0u && *written < cap) {
        for (i = 0; i < W; i++) {
            out_lo[*written * W + i] = ws_lo[i];
            out_hi[*written * W + i] = ws_hi[i];
        }
        out_unb[*written] = unb;
        (*written)++;
    }
    for (i = 0; i < W; i++) {
        last_lo[i] = ws_lo[i];
        last_hi[i] = ws_hi[i];
    }
    *last_unb = unb;
    *has_last = 1;
    (*count)++;
    return CURVE_OK;
}

CURVE_INLINE curve_status curve__box_intervals(
    int hilbert, size_t n, unsigned bits,
    const uint32_t *box_min, const uint32_t *box_max,
    uint64_t *out_lo, uint64_t *out_hi, int *out_unb,
    size_t cap, size_t *out_count, curve_range_mode mode,
    uint64_t *scratch, size_t scratch_words,
    uint32_t *axes, size_t axes_count)
{
    curve_status st;
    size_t W, key_bits, need_w, need_a, i, sp, nframes, nchild;
    size_t written = 0, count = 0;
    uint64_t *last_lo, *last_hi, *child_ws, *emit_hi, *stack;
    uint32_t *cmins, *cmaxs, *bmin, *bmax, *hscratch;
    uint32_t gmax;
    int has_last = 0, last_unb = 0, empty = 0;
    curve_cell cell, child;

    if (out_count == NULL) return CURVE_EINVAL;
    st = curve__check_dims_bits(n, bits);
    if (st != CURVE_OK) return st;
    if (n > (size_t)CURVE_RANGE_MAX_DIMS) return CURVE_EUNSUPPORTED;
    if (box_min == NULL || box_max == NULL) return CURVE_EINVAL;
    if (mode != CURVE_RANGE_EXACT && mode != CURVE_RANGE_OVERAPPROX)
        return CURVE_EINVAL;
    if (cap > 0u && (out_lo == NULL || out_hi == NULL || out_unb == NULL))
        return CURVE_EINVAL;
    if (scratch == NULL || axes == NULL) return CURVE_EINVAL;
    need_w = curve_box_intervals_scratch_words(n, bits);
    need_a = curve_box_intervals_scratch_axes(n);
    if (need_w == 0u || need_a == 0u) return CURVE_EOVERFLOW;
    if (scratch_words < need_w) return CURVE_ENOSPACE;
    if (axes_count < need_a) return CURVE_ENOSPACE;

    key_bits = n * (size_t)bits;
    W = curve__key_n_words(key_bits);
    nchild = (size_t)1 << n;
    nframes = (size_t)bits + 1u;
    last_lo = scratch;
    last_hi = scratch + W;
    child_ws = scratch + 2u * W;
    emit_hi = scratch + 3u * W;
    stack = scratch + 4u * W;
    cmins = axes;
    cmaxs = axes + n;
    bmin = axes + 2u * n;
    bmax = axes + 3u * n;
    hscratch = axes + 4u * n;

    gmax = curve__lo_mask((int)bits);
    for (i = 0; i < n; i++) {
        if (box_min[i] > box_max[i] || box_min[i] > gmax) {
            empty = 1;
            break;
        }
        bmin[i] = box_min[i];
        bmax[i] = (box_max[i] > gmax) ? gmax : box_max[i];
    }
    if (empty) {
        *out_count = 0;
        return CURVE_OK;
    }

    for (i = 0; i < W; i++) stack[i] = 0;
    stack[W] = ~0ull;
    stack[W + 1u] = 0;
    sp = 1;

    while (sp > 0u) {
        uint64_t *fr = stack + (sp - 1u) * (W + 2u);
        uint64_t *pref = fr;
        uint64_t next = fr[W];
        size_t pbits = (size_t)fr[W + 1u];
        int rel, can_merge;

        st = curve_cell_init(&cell, n, bits, pbits);
        if (st != CURVE_OK) return st;

        if (next == ~0ull) {
            if (hilbert)
                st = curve_hilbert_cell_bounds(&cell, pref, cmins, cmaxs, n,
                                               hscratch, n);
            else
                st = curve_morton_cell_bounds(&cell, pref, cmins, cmaxs, n);
            if (st != CURVE_OK) return st;
            rel = curve__aabb_relate(cmins, cmaxs, bmin, bmax, n);
            if (rel == CURVE__AABB_DISJOINT) {
                sp--;
                continue;
            }
            can_merge = 0;
            if (has_last && !last_unb) {
                st = curve_key_prefix_mask(pref, key_bits, pbits, child_ws, W);
                if (st != CURVE_OK) return st;
                if (curve_key_compare(child_ws, last_hi, W) == 0)
                    can_merge = 1;
            }
            if (rel == CURVE__AABB_CONTAINED || curve_cell_is_leaf(&cell) ||
                (mode == CURVE_RANGE_OVERAPPROX && cap > 0u &&
                 ((written >= cap && has_last) ||
                  (written + (can_merge ? 0u : 1u) >= cap)))) {
                st = curve__range_emit(pref, key_bits, pbits, W,
                                       last_lo, last_hi, &last_unb, &has_last,
                                       out_lo, out_hi, out_unb, cap,
                                       &written, &count, mode,
                                       child_ws, emit_hi);
                if (st != CURVE_OK) return st;
                sp--;
                continue;
            }
            fr[W] = 0;
            continue;
        }

        if ((size_t)next >= nchild) {
            sp--;
            continue;
        }
        st = curve_cell_child(&cell, pref, next, &child, child_ws, W);
        if (st != CURVE_OK) return st;
        fr[W] = next + 1ull;
        if (sp >= nframes) return CURVE_EOVERFLOW;
        {
            uint64_t *nf = stack + sp * (W + 2u);
            for (i = 0; i < W; i++) nf[i] = child_ws[i];
            nf[W] = ~0ull;
            nf[W + 1u] = (uint64_t)child.prefix_bits;
            sp++;
        }
    }

    *out_count = count;
    if (mode == CURVE_RANGE_EXACT && cap > 0u && count > cap)
        return CURVE_ENOSPACE;
    return CURVE_OK;
}

CURVE_INLINE curve_status curve_morton_box_intervals(
    size_t dimensions, unsigned bits_per_axis,
    const uint32_t *box_min, const uint32_t *box_max,
    uint64_t *out_lower, uint64_t *out_upper, int *out_unbounded,
    size_t interval_capacity, size_t *out_count, curve_range_mode mode,
    uint64_t *scratch, size_t scratch_words,
    uint32_t *axes_scratch, size_t axes_count)
{
    return curve__box_intervals(0, dimensions, bits_per_axis, box_min, box_max,
                                out_lower, out_upper, out_unbounded,
                                interval_capacity, out_count, mode,
                                scratch, scratch_words,
                                axes_scratch, axes_count);
}

CURVE_INLINE curve_status curve_hilbert_box_intervals(
    size_t dimensions, unsigned bits_per_axis,
    const uint32_t *box_min, const uint32_t *box_max,
    uint64_t *out_lower, uint64_t *out_upper, int *out_unbounded,
    size_t interval_capacity, size_t *out_count, curve_range_mode mode,
    uint64_t *scratch, size_t scratch_words,
    uint32_t *axes_scratch, size_t axes_count)
{
    return curve__box_intervals(1, dimensions, bits_per_axis, box_min, box_max,
                                out_lower, out_upper, out_unbounded,
                                interval_capacity, out_count, mode,
                                scratch, scratch_words,
                                axes_scratch, axes_count);
}

/* -------------------------------------------------------------------------- */
/* Pointerless prefix directory over a sorted key array                       */
/*                                                                            */
/* Binary PATRICIA tree: node.prefix_bits is the common-prefix length of      */
/* keys[begin, end). Internal nodes split at that bit. Leaves stop when the   */
/* range has `leaf_keys` or fewer keys, or the keys are identical.            */
/* nodes[0] is the root. Caller owns the node array.                          */
/* -------------------------------------------------------------------------- */

typedef struct curve_prefix_node {
    size_t prefix_bits;
    size_t begin;
    size_t end;
    size_t child0;
    size_t n_children;
} curve_prefix_node;

CURVE_INLINE size_t curve_prefix_dir_scratch_words(size_t n_keys)
{
    if (n_keys > ((size_t)-1) / 4u) return 0;
    return 4u * (n_keys + 1u);
}

CURVE_INLINE curve_status curve_prefix_dir_build(
    const uint64_t *keys, size_t n_keys, size_t stride_words,
    size_t key_bits, size_t leaf_keys,
    curve_prefix_node *nodes, size_t node_capacity, size_t *out_n_nodes,
    uint64_t *scratch, size_t scratch_words)
{
    size_t W, i, sp, used, need;
    uint64_t *stk;
    if (out_n_nodes == NULL) return CURVE_EINVAL;
    if (key_bits == 0u) return CURVE_EINVAL;
    if (n_keys > 0u && keys == NULL) return CURVE_EINVAL;
    if (n_keys > 0u && stride_words < curve__key_n_words(key_bits))
        return CURVE_EINVAL;
    if (n_keys > 0u && scratch == NULL) return CURVE_EINVAL;
    W = curve__key_n_words(key_bits);
    if (leaf_keys < 1u) leaf_keys = 1u;
    need = curve_prefix_dir_scratch_words(n_keys);
    if (need == 0u && n_keys > 0u) return CURVE_EOVERFLOW;
    if (n_keys > 0u && scratch_words < need) return CURVE_ENOSPACE;
    for (i = 1; i < n_keys; i++) {
        if (curve_key_compare(keys + (i - 1u) * stride_words,
                              keys + i * stride_words, W) > 0)
            return CURVE_EINVAL;
    }
    if (n_keys == 0u) {
        *out_n_nodes = 0;
        return CURVE_OK;
    }
    stk = scratch;
    used = 1;
    stk[0] = 0;                 /* begin */
    stk[1] = (uint64_t)n_keys;  /* end */
    stk[2] = 0;                 /* node index */
    sp = 1;
    while (sp > 0u) {
        size_t begin, end, idx, cp = 0, mid, left, right;
        const uint64_t *ka, *kb;
        sp--;
        begin = (size_t)stk[sp * 4u + 0u];
        end = (size_t)stk[sp * 4u + 1u];
        idx = (size_t)stk[sp * 4u + 2u];
        ka = keys + begin * stride_words;
        kb = keys + (end - 1u) * stride_words;
        if (curve_key_common_prefix_bits(ka, kb, key_bits, &cp) != CURVE_OK)
            return CURVE_EINVAL;
        if ((end - begin) <= leaf_keys || cp >= key_bits) {
            if (nodes != NULL && idx < node_capacity) {
                nodes[idx].prefix_bits = cp;
                nodes[idx].begin = begin;
                nodes[idx].end = end;
                nodes[idx].child0 = 0;
                nodes[idx].n_children = 0;
            }
            continue;
        }
        {
            size_t lo = begin, hi = end;
            unsigned bit;
            while (lo < hi) {
                size_t m = lo + (hi - lo) / 2u;
                bit = curve__key_get_bit(keys + m * stride_words, cp);
                if (bit == 0u) lo = m + 1u;
                else hi = m;
            }
            mid = lo;
        }
        if (mid <= begin || mid >= end) {
            /* Should not happen for a differing bit; treat as leaf. */
            if (nodes != NULL && idx < node_capacity) {
                nodes[idx].prefix_bits = cp;
                nodes[idx].begin = begin;
                nodes[idx].end = end;
                nodes[idx].child0 = 0;
                nodes[idx].n_children = 0;
            }
            continue;
        }
        left = used++;
        right = used++;
        if (nodes != NULL && idx < node_capacity) {
            nodes[idx].prefix_bits = cp;
            nodes[idx].begin = begin;
            nodes[idx].end = end;
            nodes[idx].child0 = left;
            nodes[idx].n_children = 2;
        }
        if (sp + 2u > n_keys + 1u) return CURVE_EOVERFLOW;
        stk[sp * 4u + 0u] = (uint64_t)mid;
        stk[sp * 4u + 1u] = (uint64_t)end;
        stk[sp * 4u + 2u] = (uint64_t)right;
        sp++;
        stk[sp * 4u + 0u] = (uint64_t)begin;
        stk[sp * 4u + 1u] = (uint64_t)mid;
        stk[sp * 4u + 2u] = (uint64_t)left;
        sp++;
    }
    *out_n_nodes = used;
    if (node_capacity > 0u && used > node_capacity) return CURVE_ENOSPACE;
    return CURVE_OK;
}

CURVE_INLINE curve_status curve_prefix_dir_fill_bounds(
    const curve_prefix_node *nodes, size_t n_nodes,
    const uint64_t *keys, size_t stride_words,
    size_t dimensions, unsigned bits_per_axis, int hilbert,
    uint32_t *mins, uint32_t *maxs,
    uint32_t *scratch, size_t scratch_axes)
{
    size_t i, d, k, n, W;
    curve_status st;
    uint32_t *coords;
    if (n_nodes == 0u) return CURVE_OK;
    if (nodes == NULL || keys == NULL || mins == NULL || maxs == NULL)
        return CURVE_EINVAL;
    st = curve__check_dims_bits(dimensions, bits_per_axis);
    if (st != CURVE_OK) return st;
    n = dimensions;
    W = curve_key_words_required(n, bits_per_axis);
    if (W == 0u) return CURVE_EOVERFLOW;
    if (stride_words < W) return CURVE_EINVAL;
    if (scratch == NULL || scratch_axes < n) return CURVE_ENOSPACE;
    coords = scratch;
    for (i = n_nodes; i-- > 0u; ) {
        const curve_prefix_node *nd = &nodes[i];
        uint32_t *nm = mins + i * n;
        uint32_t *nx = maxs + i * n;
        if (nd->n_children == 0u) {
            for (d = 0; d < n; d++) {
                nm[d] = ~0u;
                nx[d] = 0u;
            }
            for (k = nd->begin; k < nd->end; k++) {
                const uint64_t *kw = keys + k * stride_words;
                if (hilbert)
                    st = curve_hilbertnd_decode_words(kw, W, n, bits_per_axis,
                                                      coords, n, coords, n);
                else
                    st = curve_mortonnd_decode_words(kw, W, n, bits_per_axis,
                                                     coords, n);
                if (st != CURVE_OK) return st;
                for (d = 0; d < n; d++) {
                    if (coords[d] < nm[d]) nm[d] = coords[d];
                    if (coords[d] > nx[d]) nx[d] = coords[d];
                }
            }
            if (nd->begin >= nd->end) {
                for (d = 0; d < n; d++) nm[d] = nx[d] = 0;
            }
        } else {
            size_t c;
            for (d = 0; d < n; d++) {
                nm[d] = ~0u;
                nx[d] = 0u;
            }
            for (c = 0; c < nd->n_children; c++) {
                size_t ci = nd->child0 + c;
                const uint32_t *cm = mins + ci * n;
                const uint32_t *cx = maxs + ci * n;
                if (ci >= n_nodes) return CURVE_EINVAL;
                for (d = 0; d < n; d++) {
                    if (cm[d] < nm[d]) nm[d] = cm[d];
                    if (cx[d] > nx[d]) nx[d] = cx[d];
                }
            }
        }
    }
    return CURVE_OK;
}

CURVE_INLINE curve_status curve_prefix_dir_prefix(
    const curve_prefix_node *nodes, size_t n_nodes,
    const uint64_t *keys, size_t stride_words, size_t key_bits,
    const uint64_t *prefix, size_t prefix_bits,
    uint64_t *scratch, size_t scratch_words,
    size_t *out_begin, size_t *out_end)
{
    size_t W, ni, need;
    if (out_begin == NULL || out_end == NULL) return CURVE_EINVAL;
    if (n_nodes == 0u) {
        *out_begin = 0;
        *out_end = 0;
        return CURVE_OK;
    }
    if (nodes == NULL || keys == NULL || prefix == NULL || key_bits == 0u)
        return CURVE_EINVAL;
    if (prefix_bits > key_bits) return CURVE_EINVAL;
    W = curve__key_n_words(key_bits);
    if (stride_words < W) return CURVE_EINVAL;
    need = 2u * W;
    if (scratch == NULL || scratch_words < need) return CURVE_ENOSPACE;
    ni = 0;
    while (nodes[ni].n_children > 0u) {
        size_t m, c, next;
        int match = 0;
        unsigned qbit, cbit;
        m = (prefix_bits < nodes[ni].prefix_bits) ? prefix_bits
                                                  : nodes[ni].prefix_bits;
        if (curve_key_prefix_matches(keys + nodes[ni].begin * stride_words,
                                     prefix, key_bits, m, &match) != CURVE_OK)
            return CURVE_EINVAL;
        if (!match) {
            *out_begin = 0;
            *out_end = 0;
            return CURVE_OK;
        }
        if (prefix_bits <= nodes[ni].prefix_bits) break;
        qbit = curve__key_get_bit(prefix, nodes[ni].prefix_bits);
        next = n_nodes;
        for (c = 0; c < nodes[ni].n_children; c++) {
            size_t ci = nodes[ni].child0 + c;
            if (ci >= n_nodes) return CURVE_EINVAL;
            cbit = curve__key_get_bit(keys + nodes[ci].begin * stride_words,
                                      nodes[ni].prefix_bits);
            if (cbit == qbit) {
                next = ci;
                break;
            }
        }
        if (next >= n_nodes) {
            *out_begin = 0;
            *out_end = 0;
            return CURVE_OK;
        }
        ni = next;
    }
    {
        curve_status st;
        size_t b, e, base;
        base = nodes[ni].begin;
        st = curve_key_prefix_range(keys + base * stride_words,
                                    nodes[ni].end - base, stride_words,
                                    prefix, key_bits, prefix_bits,
                                    scratch, scratch_words, &b, &e);
        if (st != CURVE_OK) return st;
        *out_begin = base + b;
        *out_end = base + e;
        return CURVE_OK;
    }
}

CURVE_INLINE curve_status curve_prefix_dir_box(
    const curve_prefix_node *nodes, size_t n_nodes,
    const uint32_t *node_mins, const uint32_t *node_maxs, size_t dimensions,
    const uint32_t *box_min, const uint32_t *box_max,
    const uint64_t *keys, size_t stride_words, unsigned bits_per_axis,
    int hilbert, uint32_t *axes_scratch, size_t axes_count,
    uint64_t *scratch, size_t scratch_words,
    size_t *out_begin, size_t *out_end, size_t range_capacity,
    size_t *out_nranges)
{
    size_t sp, n, W, written = 0, count = 0, last_e = 0;
    int has_r = 0;
    uint64_t *stk;
    uint32_t *coords;
    if (out_nranges == NULL) return CURVE_EINVAL;
    if (n_nodes == 0u) {
        *out_nranges = 0;
        return CURVE_OK;
    }
    if (nodes == NULL || node_mins == NULL || node_maxs == NULL)
        return CURVE_EINVAL;
    if (box_min == NULL || box_max == NULL) return CURVE_EINVAL;
    n = dimensions;
    if (n < 1u) return CURVE_EINVAL;
    if (range_capacity > 0u && (out_begin == NULL || out_end == NULL))
        return CURVE_EINVAL;
    if (scratch == NULL || scratch_words < n_nodes) return CURVE_ENOSPACE;
    W = 0;
    coords = NULL;
    if (keys != NULL) {
        curve_status st = curve__check_dims_bits(n, bits_per_axis);
        if (st != CURVE_OK) return st;
        W = curve_key_words_required(n, bits_per_axis);
        if (stride_words < W) return CURVE_EINVAL;
        if (axes_scratch == NULL || axes_count < n) return CURVE_ENOSPACE;
        coords = axes_scratch;
    }
    stk = scratch;
    stk[0] = 0;
    sp = 1;
    while (sp > 0u) {
        size_t ni, rel, c;
        const uint32_t *nm, *nx;
        sp--;
        ni = (size_t)stk[sp];
        if (ni >= n_nodes) return CURVE_EINVAL;
        nm = node_mins + ni * n;
        nx = node_maxs + ni * n;
        rel = (size_t)curve__aabb_relate(nm, nx, box_min, box_max, n);
        if (rel == (size_t)CURVE__AABB_DISJOINT) continue;
        if (rel == (size_t)CURVE__AABB_CONTAINED ||
            (nodes[ni].n_children == 0u && keys == NULL)) {
            size_t b = nodes[ni].begin, e = nodes[ni].end;
            if (has_r && last_e == b) {
                last_e = e;
                if (range_capacity > 0u && written > 0u)
                    out_end[written - 1u] = e;
            } else {
                if (range_capacity > 0u && written < range_capacity) {
                    out_begin[written] = b;
                    out_end[written] = e;
                    written++;
                }
                last_e = e;
                has_r = 1;
                count++;
            }
            continue;
        }
        if (nodes[ni].n_children == 0u && keys != NULL) {
            size_t k;
            for (k = nodes[ni].begin; k < nodes[ni].end; k++) {
                const uint64_t *kw = keys + k * stride_words;
                int inside = 1;
                curve_status st;
                size_t d;
                if (hilbert)
                    st = curve_hilbertnd_decode_words(kw, W, n, bits_per_axis,
                                                      coords, n, coords, n);
                else
                    st = curve_mortonnd_decode_words(kw, W, n, bits_per_axis,
                                                     coords, n);
                if (st != CURVE_OK) return st;
                for (d = 0; d < n; d++) {
                    if (coords[d] < box_min[d] || coords[d] > box_max[d]) {
                        inside = 0;
                        break;
                    }
                }
                if (!inside) continue;
                if (has_r && last_e == k) {
                    last_e = k + 1u;
                    if (range_capacity > 0u && written > 0u)
                        out_end[written - 1u] = last_e;
                } else {
                    if (range_capacity > 0u && written < range_capacity) {
                        out_begin[written] = k;
                        out_end[written] = k + 1u;
                        written++;
                    }
                    last_e = k + 1u;
                    has_r = 1;
                    count++;
                }
            }
            continue;
        }
        for (c = nodes[ni].n_children; c-- > 0u; ) {
            if (sp >= n_nodes) return CURVE_EOVERFLOW;
            stk[sp++] = (uint64_t)(nodes[ni].child0 + c);
        }
    }
    *out_nranges = count;
    if (range_capacity > 0u && count > range_capacity)
        return CURVE_ENOSPACE;
    return CURVE_OK;
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
/* present is used; a scalar tail handles leftovers.  32-bit Morton 2D and    */
/* 32-bit Hilbert 2D use NEON / AVX2 / AVX-512.  64-bit Morton 2D/3D and      */
/* 64-bit Hilbert 2D batches use AVX2 / AVX-512 only (2-wide NEON SWAR loses  */
/* to the compiler's scalar loop on aarch64, and Hilbert-64 d2xy is worse).   */
/* n-D Hilbert / Skilling stays scalar: the transform couples axes with       */
/* data-dependent swaps, and 256-D cost is dominated by that walk.  Scalar    */
/* BMI2 PDEP/PEXT still applies to the single-element APIs above.             */
/*                                                                            */
/* curve_compiled_isa() names the compile-time path.  curve_simd_isa() is the */
/* historical alias.  Runtime selection uses a caller-owned curve_dispatch    */
/* (no process globals).  CURVE_RUNTIME_DISPATCH compiles extra x86 kernels.  */
/* -------------------------------------------------------------------------- */

CURVE_INLINE const char *curve_compiled_isa(void)
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

CURVE_INLINE const char *curve_simd_isa(void)
{
    return curve_compiled_isa();
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

#if defined(CURVE_HAVE_AVX2) || defined(CURVE_DISPATCH_X86)
CURVE_ATTR_AVX2 CURVE_INLINE __m256i curve__part1by1_avx2(__m256i x)
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
CURVE_ATTR_AVX2 CURVE_INLINE __m256i curve__compact1by1_avx2(__m256i x)
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
CURVE_ATTR_AVX2 CURVE_INLINE __m256i curve__part1by1_avx2_64(__m256i x)
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
CURVE_ATTR_AVX2 CURVE_INLINE __m256i curve__compact1by1_avx2_64(__m256i x)
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
CURVE_ATTR_AVX2 CURVE_INLINE __m256i curve__part1by2_avx2_64(__m256i x)
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
CURVE_ATTR_AVX2 CURVE_INLINE __m256i curve__compact1by2_avx2_64(__m256i x)
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
CURVE_ATTR_AVX2 CURVE_INLINE __m128i curve__pack64lo32_avx2(__m256i v)
{
    /* Low 32 bits of each 64-bit lane -> 4x uint32. */
    v = _mm256_shuffle_epi32(v, _MM_SHUFFLE(2, 0, 2, 0));
    return _mm_unpacklo_epi64(_mm256_castsi256_si128(v),
                              _mm256_extracti128_si256(v, 1));
}
#endif /* CURVE_HAVE_AVX2 */

#if defined(CURVE_HAVE_AVX512) || defined(CURVE_DISPATCH_X86)
CURVE_ATTR_AVX512 CURVE_INLINE __m512i curve__part1by1_avx512(__m512i x)
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
CURVE_ATTR_AVX512 CURVE_INLINE __m512i curve__compact1by1_avx512(__m512i x)
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
CURVE_ATTR_AVX512 CURVE_INLINE __m512i curve__part1by1_avx512_64(__m512i x)
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
CURVE_ATTR_AVX512 CURVE_INLINE __m512i curve__compact1by1_avx512_64(__m512i x)
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
CURVE_ATTR_AVX512 CURVE_INLINE __m512i curve__part1by2_avx512_64(__m512i x)
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
CURVE_ATTR_AVX512 CURVE_INLINE __m512i curve__compact1by2_avx512_64(__m512i x)
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

/* Hilbert 2D SIMD cores — same bitwise map as curve_hilbert_xy2d / d2xy and
 * the 64-bit variants.  `bits` is already clamped by the batch wrapper. */

#if defined(CURVE_HAVE_NEON)
CURVE_INLINE uint32x4_t curve__prefixscan_u32x4(uint32x4_t x)
{
    x = veorq_u32(vshrq_n_u32(x, 8), x);
    x = veorq_u32(vshrq_n_u32(x, 4), x);
    x = veorq_u32(vshrq_n_u32(x, 2), x);
    x = veorq_u32(vshrq_n_u32(x, 1), x);
    return x;
}
CURVE_INLINE uint32x4_t curve__hilbert_xy2d_u32x4(uint32x4_t x, uint32x4_t y, int bits)
{
    const uint32x4_t F = vdupq_n_u32(0xFFFFu);
    uint32x4_t a, b, c, d, A, B, C, D, i0, i1;
    x = vandq_u32(x, vdupq_n_u32(curve__lo_mask(bits)));
    y = vandq_u32(y, vdupq_n_u32(curve__lo_mask(bits)));
    x = vshlq_u32(x, vdupq_n_s32(16 - bits));
    y = vshlq_u32(y, vdupq_n_s32(16 - bits));
    a = veorq_u32(x, y);
    b = veorq_u32(F, a);
    c = veorq_u32(F, vorrq_u32(x, y));
    d = vandq_u32(x, veorq_u32(y, F));
    A = vorrq_u32(a, vshrq_n_u32(b, 1));
    B = veorq_u32(vshrq_n_u32(a, 1), a);
    C = veorq_u32(veorq_u32(vshrq_n_u32(c, 1), vandq_u32(b, vshrq_n_u32(d, 1))), c);
    D = veorq_u32(veorq_u32(vandq_u32(a, vshrq_n_u32(c, 1)), vshrq_n_u32(d, 1)), d);
    a = A; b = B; c = C; d = D;
    A = veorq_u32(vandq_u32(a, vshrq_n_u32(a, 2)), vandq_u32(b, vshrq_n_u32(b, 2)));
    B = veorq_u32(vandq_u32(a, vshrq_n_u32(b, 2)), vandq_u32(b, vshrq_n_u32(veorq_u32(a, b), 2)));
    C = veorq_u32(C, veorq_u32(vandq_u32(a, vshrq_n_u32(c, 2)), vandq_u32(b, vshrq_n_u32(d, 2))));
    D = veorq_u32(D, veorq_u32(vandq_u32(b, vshrq_n_u32(c, 2)), vandq_u32(veorq_u32(a, b), vshrq_n_u32(d, 2))));
    a = A; b = B; c = C; d = D;
    A = veorq_u32(vandq_u32(a, vshrq_n_u32(a, 4)), vandq_u32(b, vshrq_n_u32(b, 4)));
    B = veorq_u32(vandq_u32(a, vshrq_n_u32(b, 4)), vandq_u32(b, vshrq_n_u32(veorq_u32(a, b), 4)));
    C = veorq_u32(C, veorq_u32(vandq_u32(a, vshrq_n_u32(c, 4)), vandq_u32(b, vshrq_n_u32(d, 4))));
    D = veorq_u32(D, veorq_u32(vandq_u32(b, vshrq_n_u32(c, 4)), vandq_u32(veorq_u32(a, b), vshrq_n_u32(d, 4))));
    a = A; b = B; c = C; d = D;
    C = veorq_u32(C, veorq_u32(vandq_u32(a, vshrq_n_u32(c, 8)), vandq_u32(b, vshrq_n_u32(d, 8))));
    D = veorq_u32(D, veorq_u32(vandq_u32(b, vshrq_n_u32(c, 8)), vandq_u32(veorq_u32(a, b), vshrq_n_u32(d, 8))));
    a = veorq_u32(C, vshrq_n_u32(C, 1));
    b = veorq_u32(D, vshrq_n_u32(D, 1));
    i0 = veorq_u32(x, y);
    i1 = vorrq_u32(b, veorq_u32(F, vorrq_u32(i0, a)));
    return vshlq_u32(vorrq_u32(vshlq_n_u32(curve__part1by1_u32x4(i1), 1),
                               curve__part1by1_u32x4(i0)),
                     vdupq_n_s32(-(32 - 2 * bits)));
}
CURVE_INLINE void curve__hilbert_d2xy_u32x4(uint32x4_t di, int bits,
                                           uint32x4_t *CURVE_RESTRICT x,
                                           uint32x4_t *CURVE_RESTRICT y)
{
    const uint32x4_t F = vdupq_n_u32(0xFFFFu);
    uint32x4_t i, i0, i1, t0, t1, a;
    i = vshlq_u32(di, vdupq_n_s32(32 - 2 * bits));
    i0 = curve__compact1by1_u32x4(i);
    i1 = curve__compact1by1_u32x4(vshrq_n_u32(i, 1));
    t0 = veorq_u32(vorrq_u32(i0, i1), F);
    t1 = vandq_u32(i0, i1);
    a = vorrq_u32(vandq_u32(veorq_u32(i0, F), curve__prefixscan_u32x4(t1)),
                  vandq_u32(i0, curve__prefixscan_u32x4(t0)));
    *x = vshlq_u32(veorq_u32(a, i1), vdupq_n_s32(-(16 - bits)));
    *y = vshlq_u32(veorq_u32(veorq_u32(a, i0), i1), vdupq_n_s32(-(16 - bits)));
}
#endif /* CURVE_HAVE_NEON */

#if defined(CURVE_HAVE_AVX2) || defined(CURVE_DISPATCH_X86)
CURVE_ATTR_AVX2 CURVE_INLINE __m256i curve__prefixscan_avx2(__m256i x)
{
    x = _mm256_xor_si256(_mm256_srli_epi32(x, 8), x);
    x = _mm256_xor_si256(_mm256_srli_epi32(x, 4), x);
    x = _mm256_xor_si256(_mm256_srli_epi32(x, 2), x);
    x = _mm256_xor_si256(_mm256_srli_epi32(x, 1), x);
    return x;
}
CURVE_ATTR_AVX2 CURVE_INLINE __m256i curve__prefixscan_avx2_64(__m256i x)
{
    x = _mm256_xor_si256(_mm256_srli_epi64(x, 16), x);
    x = _mm256_xor_si256(_mm256_srli_epi64(x, 8), x);
    x = _mm256_xor_si256(_mm256_srli_epi64(x, 4), x);
    x = _mm256_xor_si256(_mm256_srli_epi64(x, 2), x);
    x = _mm256_xor_si256(_mm256_srli_epi64(x, 1), x);
    return x;
}
CURVE_ATTR_AVX2 CURVE_INLINE __m256i curve__hilbert_xy2d_avx2(__m256i x, __m256i y, int bits)
{
    const __m256i F = _mm256_set1_epi32((int)0xFFFFu);
    __m256i a, b, c, d, A, B, C, D, i0, i1;
    x = _mm256_and_si256(x, _mm256_set1_epi32((int)curve__lo_mask(bits)));
    y = _mm256_and_si256(y, _mm256_set1_epi32((int)curve__lo_mask(bits)));
    x = _mm256_sllv_epi32(x, _mm256_set1_epi32(16 - bits));
    y = _mm256_sllv_epi32(y, _mm256_set1_epi32(16 - bits));
    a = _mm256_xor_si256(x, y);
    b = _mm256_xor_si256(F, a);
    c = _mm256_xor_si256(F, _mm256_or_si256(x, y));
    d = _mm256_and_si256(x, _mm256_xor_si256(y, F));
    A = _mm256_or_si256(a, _mm256_srli_epi32(b, 1));
    B = _mm256_xor_si256(_mm256_srli_epi32(a, 1), a);
    C = _mm256_xor_si256(_mm256_xor_si256(_mm256_srli_epi32(c, 1),
                          _mm256_and_si256(b, _mm256_srli_epi32(d, 1))), c);
    D = _mm256_xor_si256(_mm256_xor_si256(_mm256_and_si256(a, _mm256_srli_epi32(c, 1)),
                          _mm256_srli_epi32(d, 1)), d);
    a = A; b = B; c = C; d = D;
    A = _mm256_xor_si256(_mm256_and_si256(a, _mm256_srli_epi32(a, 2)),
                         _mm256_and_si256(b, _mm256_srli_epi32(b, 2)));
    B = _mm256_xor_si256(_mm256_and_si256(a, _mm256_srli_epi32(b, 2)),
                         _mm256_and_si256(b, _mm256_srli_epi32(_mm256_xor_si256(a, b), 2)));
    C = _mm256_xor_si256(C, _mm256_xor_si256(_mm256_and_si256(a, _mm256_srli_epi32(c, 2)),
                         _mm256_and_si256(b, _mm256_srli_epi32(d, 2))));
    D = _mm256_xor_si256(D, _mm256_xor_si256(_mm256_and_si256(b, _mm256_srli_epi32(c, 2)),
                         _mm256_and_si256(_mm256_xor_si256(a, b), _mm256_srli_epi32(d, 2))));
    a = A; b = B; c = C; d = D;
    A = _mm256_xor_si256(_mm256_and_si256(a, _mm256_srli_epi32(a, 4)),
                         _mm256_and_si256(b, _mm256_srli_epi32(b, 4)));
    B = _mm256_xor_si256(_mm256_and_si256(a, _mm256_srli_epi32(b, 4)),
                         _mm256_and_si256(b, _mm256_srli_epi32(_mm256_xor_si256(a, b), 4)));
    C = _mm256_xor_si256(C, _mm256_xor_si256(_mm256_and_si256(a, _mm256_srli_epi32(c, 4)),
                         _mm256_and_si256(b, _mm256_srli_epi32(d, 4))));
    D = _mm256_xor_si256(D, _mm256_xor_si256(_mm256_and_si256(b, _mm256_srli_epi32(c, 4)),
                         _mm256_and_si256(_mm256_xor_si256(a, b), _mm256_srli_epi32(d, 4))));
    a = A; b = B; c = C; d = D;
    C = _mm256_xor_si256(C, _mm256_xor_si256(_mm256_and_si256(a, _mm256_srli_epi32(c, 8)),
                         _mm256_and_si256(b, _mm256_srli_epi32(d, 8))));
    D = _mm256_xor_si256(D, _mm256_xor_si256(_mm256_and_si256(b, _mm256_srli_epi32(c, 8)),
                         _mm256_and_si256(_mm256_xor_si256(a, b), _mm256_srli_epi32(d, 8))));
    a = _mm256_xor_si256(C, _mm256_srli_epi32(C, 1));
    b = _mm256_xor_si256(D, _mm256_srli_epi32(D, 1));
    i0 = _mm256_xor_si256(x, y);
    i1 = _mm256_or_si256(b, _mm256_xor_si256(F, _mm256_or_si256(i0, a)));
    return _mm256_srlv_epi32(_mm256_or_si256(_mm256_slli_epi32(curve__part1by1_avx2(i1), 1),
                                             curve__part1by1_avx2(i0)),
                             _mm256_set1_epi32(32 - 2 * bits));
}
CURVE_ATTR_AVX2 CURVE_INLINE void curve__hilbert_d2xy_avx2(__m256i di, int bits,
                                          __m256i *CURVE_RESTRICT x,
                                          __m256i *CURVE_RESTRICT y)
{
    const __m256i F = _mm256_set1_epi32((int)0xFFFFu);
    __m256i i, i0, i1, t0, t1, a;
    i = _mm256_sllv_epi32(di, _mm256_set1_epi32(32 - 2 * bits));
    i0 = curve__compact1by1_avx2(i);
    i1 = curve__compact1by1_avx2(_mm256_srli_epi32(i, 1));
    t0 = _mm256_xor_si256(_mm256_or_si256(i0, i1), F);
    t1 = _mm256_and_si256(i0, i1);
    a = _mm256_or_si256(_mm256_and_si256(_mm256_xor_si256(i0, F), curve__prefixscan_avx2(t1)),
                        _mm256_and_si256(i0, curve__prefixscan_avx2(t0)));
    *x = _mm256_srlv_epi32(_mm256_xor_si256(a, i1), _mm256_set1_epi32(16 - bits));
    *y = _mm256_srlv_epi32(_mm256_xor_si256(_mm256_xor_si256(a, i0), i1),
                           _mm256_set1_epi32(16 - bits));
}
CURVE_ATTR_AVX2 CURVE_INLINE __m256i curve__hilbert2_64_xy2d_avx2(__m256i X, __m256i Y, int bits)
{
    const __m256i M = _mm256_set1_epi64x((long long)0xFFFFFFFFull);
    __m256i a, b, c, d, A, B, C, D, i0, i1;
    X = _mm256_and_si256(X, _mm256_set1_epi64x((long long)curve__lo_mask(bits)));
    Y = _mm256_and_si256(Y, _mm256_set1_epi64x((long long)curve__lo_mask(bits)));
    X = _mm256_sllv_epi64(X, _mm256_set1_epi64x((long long)(32 - bits)));
    Y = _mm256_sllv_epi64(Y, _mm256_set1_epi64x((long long)(32 - bits)));
    a = _mm256_xor_si256(X, Y);
    b = _mm256_xor_si256(M, a);
    c = _mm256_xor_si256(M, _mm256_or_si256(X, Y));
    d = _mm256_and_si256(X, _mm256_xor_si256(Y, M));
    A = _mm256_or_si256(a, _mm256_srli_epi64(b, 1));
    B = _mm256_xor_si256(_mm256_srli_epi64(a, 1), a);
    C = _mm256_xor_si256(_mm256_xor_si256(_mm256_srli_epi64(c, 1),
                          _mm256_and_si256(b, _mm256_srli_epi64(d, 1))), c);
    D = _mm256_xor_si256(_mm256_xor_si256(_mm256_and_si256(a, _mm256_srli_epi64(c, 1)),
                          _mm256_srli_epi64(d, 1)), d);
    a = A; b = B; c = C; d = D;
    A = _mm256_xor_si256(_mm256_and_si256(a, _mm256_srli_epi64(a, 2)),
                         _mm256_and_si256(b, _mm256_srli_epi64(b, 2)));
    B = _mm256_xor_si256(_mm256_and_si256(a, _mm256_srli_epi64(b, 2)),
                         _mm256_and_si256(b, _mm256_srli_epi64(_mm256_xor_si256(a, b), 2)));
    C = _mm256_xor_si256(C, _mm256_xor_si256(_mm256_and_si256(a, _mm256_srli_epi64(c, 2)),
                         _mm256_and_si256(b, _mm256_srli_epi64(d, 2))));
    D = _mm256_xor_si256(D, _mm256_xor_si256(_mm256_and_si256(b, _mm256_srli_epi64(c, 2)),
                         _mm256_and_si256(_mm256_xor_si256(a, b), _mm256_srli_epi64(d, 2))));
    a = A; b = B; c = C; d = D;
    A = _mm256_xor_si256(_mm256_and_si256(a, _mm256_srli_epi64(a, 4)),
                         _mm256_and_si256(b, _mm256_srli_epi64(b, 4)));
    B = _mm256_xor_si256(_mm256_and_si256(a, _mm256_srli_epi64(b, 4)),
                         _mm256_and_si256(b, _mm256_srli_epi64(_mm256_xor_si256(a, b), 4)));
    C = _mm256_xor_si256(C, _mm256_xor_si256(_mm256_and_si256(a, _mm256_srli_epi64(c, 4)),
                         _mm256_and_si256(b, _mm256_srli_epi64(d, 4))));
    D = _mm256_xor_si256(D, _mm256_xor_si256(_mm256_and_si256(b, _mm256_srli_epi64(c, 4)),
                         _mm256_and_si256(_mm256_xor_si256(a, b), _mm256_srli_epi64(d, 4))));
    a = A; b = B; c = C; d = D;
    A = _mm256_xor_si256(_mm256_and_si256(a, _mm256_srli_epi64(a, 8)),
                         _mm256_and_si256(b, _mm256_srli_epi64(b, 8)));
    B = _mm256_xor_si256(_mm256_and_si256(a, _mm256_srli_epi64(b, 8)),
                         _mm256_and_si256(b, _mm256_srli_epi64(_mm256_xor_si256(a, b), 8)));
    C = _mm256_xor_si256(C, _mm256_xor_si256(_mm256_and_si256(a, _mm256_srli_epi64(c, 8)),
                         _mm256_and_si256(b, _mm256_srli_epi64(d, 8))));
    D = _mm256_xor_si256(D, _mm256_xor_si256(_mm256_and_si256(b, _mm256_srli_epi64(c, 8)),
                         _mm256_and_si256(_mm256_xor_si256(a, b), _mm256_srli_epi64(d, 8))));
    a = A; b = B; c = C; d = D;
    C = _mm256_xor_si256(C, _mm256_xor_si256(_mm256_and_si256(a, _mm256_srli_epi64(c, 16)),
                         _mm256_and_si256(b, _mm256_srli_epi64(d, 16))));
    D = _mm256_xor_si256(D, _mm256_xor_si256(_mm256_and_si256(b, _mm256_srli_epi64(c, 16)),
                         _mm256_and_si256(_mm256_xor_si256(a, b), _mm256_srli_epi64(d, 16))));
    a = _mm256_xor_si256(C, _mm256_srli_epi64(C, 1));
    b = _mm256_xor_si256(D, _mm256_srli_epi64(D, 1));
    i0 = _mm256_xor_si256(X, Y);
    i1 = _mm256_or_si256(b, _mm256_xor_si256(M, _mm256_or_si256(i0, a)));
    return _mm256_srlv_epi64(_mm256_or_si256(_mm256_slli_epi64(curve__part1by1_avx2_64(i1), 1),
                                             curve__part1by1_avx2_64(i0)),
                             _mm256_set1_epi64x((long long)(64 - 2 * bits)));
}
CURVE_ATTR_AVX2 CURVE_INLINE void curve__hilbert2_64_d2xy_avx2(__m256i di, int bits,
                                              __m256i *CURVE_RESTRICT x,
                                              __m256i *CURVE_RESTRICT y)
{
    const __m256i M = _mm256_set1_epi64x((long long)0xFFFFFFFFull);
    __m256i i, i0, i1, t0, t1, a;
    i = _mm256_sllv_epi64(di, _mm256_set1_epi64x((long long)(64 - 2 * bits)));
    i0 = curve__compact1by1_avx2_64(i);
    i1 = curve__compact1by1_avx2_64(_mm256_srli_epi64(i, 1));
    t0 = _mm256_xor_si256(_mm256_or_si256(i0, i1), M);
    t1 = _mm256_and_si256(i0, i1);
    a = _mm256_or_si256(_mm256_and_si256(_mm256_xor_si256(i0, M), curve__prefixscan_avx2_64(t1)),
                        _mm256_and_si256(i0, curve__prefixscan_avx2_64(t0)));
    *x = _mm256_srlv_epi64(_mm256_xor_si256(a, i1),
                           _mm256_set1_epi64x((long long)(32 - bits)));
    *y = _mm256_srlv_epi64(_mm256_xor_si256(_mm256_xor_si256(a, i0), i1),
                           _mm256_set1_epi64x((long long)(32 - bits)));
}
#endif /* CURVE_HAVE_AVX2 */

#if defined(CURVE_HAVE_AVX512) || defined(CURVE_DISPATCH_X86)
CURVE_ATTR_AVX512 CURVE_INLINE __m512i curve__prefixscan_avx512(__m512i x)
{
    x = _mm512_xor_si512(_mm512_srli_epi32(x, 8), x);
    x = _mm512_xor_si512(_mm512_srli_epi32(x, 4), x);
    x = _mm512_xor_si512(_mm512_srli_epi32(x, 2), x);
    x = _mm512_xor_si512(_mm512_srli_epi32(x, 1), x);
    return x;
}
CURVE_ATTR_AVX512 CURVE_INLINE __m512i curve__prefixscan_avx512_64(__m512i x)
{
    x = _mm512_xor_si512(_mm512_srli_epi64(x, 16), x);
    x = _mm512_xor_si512(_mm512_srli_epi64(x, 8), x);
    x = _mm512_xor_si512(_mm512_srli_epi64(x, 4), x);
    x = _mm512_xor_si512(_mm512_srli_epi64(x, 2), x);
    x = _mm512_xor_si512(_mm512_srli_epi64(x, 1), x);
    return x;
}
CURVE_ATTR_AVX512 CURVE_INLINE __m512i curve__hilbert_xy2d_avx512(__m512i x, __m512i y, int bits)
{
    const __m512i F = _mm512_set1_epi32((int)0xFFFFu);
    __m512i a, b, c, d, A, B, C, D, i0, i1;
    x = _mm512_and_si512(x, _mm512_set1_epi32((int)curve__lo_mask(bits)));
    y = _mm512_and_si512(y, _mm512_set1_epi32((int)curve__lo_mask(bits)));
    x = _mm512_sllv_epi32(x, _mm512_set1_epi32(16 - bits));
    y = _mm512_sllv_epi32(y, _mm512_set1_epi32(16 - bits));
    a = _mm512_xor_si512(x, y);
    b = _mm512_xor_si512(F, a);
    c = _mm512_xor_si512(F, _mm512_or_si512(x, y));
    d = _mm512_and_si512(x, _mm512_xor_si512(y, F));
    A = _mm512_or_si512(a, _mm512_srli_epi32(b, 1));
    B = _mm512_xor_si512(_mm512_srli_epi32(a, 1), a);
    C = _mm512_xor_si512(_mm512_xor_si512(_mm512_srli_epi32(c, 1),
                          _mm512_and_si512(b, _mm512_srli_epi32(d, 1))), c);
    D = _mm512_xor_si512(_mm512_xor_si512(_mm512_and_si512(a, _mm512_srli_epi32(c, 1)),
                          _mm512_srli_epi32(d, 1)), d);
    a = A; b = B; c = C; d = D;
    A = _mm512_xor_si512(_mm512_and_si512(a, _mm512_srli_epi32(a, 2)),
                         _mm512_and_si512(b, _mm512_srli_epi32(b, 2)));
    B = _mm512_xor_si512(_mm512_and_si512(a, _mm512_srli_epi32(b, 2)),
                         _mm512_and_si512(b, _mm512_srli_epi32(_mm512_xor_si512(a, b), 2)));
    C = _mm512_xor_si512(C, _mm512_xor_si512(_mm512_and_si512(a, _mm512_srli_epi32(c, 2)),
                         _mm512_and_si512(b, _mm512_srli_epi32(d, 2))));
    D = _mm512_xor_si512(D, _mm512_xor_si512(_mm512_and_si512(b, _mm512_srli_epi32(c, 2)),
                         _mm512_and_si512(_mm512_xor_si512(a, b), _mm512_srli_epi32(d, 2))));
    a = A; b = B; c = C; d = D;
    A = _mm512_xor_si512(_mm512_and_si512(a, _mm512_srli_epi32(a, 4)),
                         _mm512_and_si512(b, _mm512_srli_epi32(b, 4)));
    B = _mm512_xor_si512(_mm512_and_si512(a, _mm512_srli_epi32(b, 4)),
                         _mm512_and_si512(b, _mm512_srli_epi32(_mm512_xor_si512(a, b), 4)));
    C = _mm512_xor_si512(C, _mm512_xor_si512(_mm512_and_si512(a, _mm512_srli_epi32(c, 4)),
                         _mm512_and_si512(b, _mm512_srli_epi32(d, 4))));
    D = _mm512_xor_si512(D, _mm512_xor_si512(_mm512_and_si512(b, _mm512_srli_epi32(c, 4)),
                         _mm512_and_si512(_mm512_xor_si512(a, b), _mm512_srli_epi32(d, 4))));
    a = A; b = B; c = C; d = D;
    C = _mm512_xor_si512(C, _mm512_xor_si512(_mm512_and_si512(a, _mm512_srli_epi32(c, 8)),
                         _mm512_and_si512(b, _mm512_srli_epi32(d, 8))));
    D = _mm512_xor_si512(D, _mm512_xor_si512(_mm512_and_si512(b, _mm512_srli_epi32(c, 8)),
                         _mm512_and_si512(_mm512_xor_si512(a, b), _mm512_srli_epi32(d, 8))));
    a = _mm512_xor_si512(C, _mm512_srli_epi32(C, 1));
    b = _mm512_xor_si512(D, _mm512_srli_epi32(D, 1));
    i0 = _mm512_xor_si512(x, y);
    i1 = _mm512_or_si512(b, _mm512_xor_si512(F, _mm512_or_si512(i0, a)));
    return _mm512_srlv_epi32(_mm512_or_si512(_mm512_slli_epi32(curve__part1by1_avx512(i1), 1),
                                             curve__part1by1_avx512(i0)),
                             _mm512_set1_epi32(32 - 2 * bits));
}
CURVE_ATTR_AVX512 CURVE_INLINE void curve__hilbert_d2xy_avx512(__m512i di, int bits,
                                            __m512i *CURVE_RESTRICT x,
                                            __m512i *CURVE_RESTRICT y)
{
    const __m512i F = _mm512_set1_epi32((int)0xFFFFu);
    __m512i i, i0, i1, t0, t1, a;
    i = _mm512_sllv_epi32(di, _mm512_set1_epi32(32 - 2 * bits));
    i0 = curve__compact1by1_avx512(i);
    i1 = curve__compact1by1_avx512(_mm512_srli_epi32(i, 1));
    t0 = _mm512_xor_si512(_mm512_or_si512(i0, i1), F);
    t1 = _mm512_and_si512(i0, i1);
    a = _mm512_or_si512(_mm512_and_si512(_mm512_xor_si512(i0, F), curve__prefixscan_avx512(t1)),
                        _mm512_and_si512(i0, curve__prefixscan_avx512(t0)));
    *x = _mm512_srlv_epi32(_mm512_xor_si512(a, i1), _mm512_set1_epi32(16 - bits));
    *y = _mm512_srlv_epi32(_mm512_xor_si512(_mm512_xor_si512(a, i0), i1),
                           _mm512_set1_epi32(16 - bits));
}
CURVE_ATTR_AVX512 CURVE_INLINE __m512i curve__hilbert2_64_xy2d_avx512(__m512i X, __m512i Y, int bits)
{
    const __m512i M = _mm512_set1_epi64((long long)0xFFFFFFFFull);
    __m512i a, b, c, d, A, B, C, D, i0, i1;
    X = _mm512_and_si512(X, _mm512_set1_epi64((long long)curve__lo_mask(bits)));
    Y = _mm512_and_si512(Y, _mm512_set1_epi64((long long)curve__lo_mask(bits)));
    X = _mm512_sllv_epi64(X, _mm512_set1_epi64((long long)(32 - bits)));
    Y = _mm512_sllv_epi64(Y, _mm512_set1_epi64((long long)(32 - bits)));
    a = _mm512_xor_si512(X, Y);
    b = _mm512_xor_si512(M, a);
    c = _mm512_xor_si512(M, _mm512_or_si512(X, Y));
    d = _mm512_and_si512(X, _mm512_xor_si512(Y, M));
    A = _mm512_or_si512(a, _mm512_srli_epi64(b, 1));
    B = _mm512_xor_si512(_mm512_srli_epi64(a, 1), a);
    C = _mm512_xor_si512(_mm512_xor_si512(_mm512_srli_epi64(c, 1),
                          _mm512_and_si512(b, _mm512_srli_epi64(d, 1))), c);
    D = _mm512_xor_si512(_mm512_xor_si512(_mm512_and_si512(a, _mm512_srli_epi64(c, 1)),
                          _mm512_srli_epi64(d, 1)), d);
    a = A; b = B; c = C; d = D;
    A = _mm512_xor_si512(_mm512_and_si512(a, _mm512_srli_epi64(a, 2)),
                         _mm512_and_si512(b, _mm512_srli_epi64(b, 2)));
    B = _mm512_xor_si512(_mm512_and_si512(a, _mm512_srli_epi64(b, 2)),
                         _mm512_and_si512(b, _mm512_srli_epi64(_mm512_xor_si512(a, b), 2)));
    C = _mm512_xor_si512(C, _mm512_xor_si512(_mm512_and_si512(a, _mm512_srli_epi64(c, 2)),
                         _mm512_and_si512(b, _mm512_srli_epi64(d, 2))));
    D = _mm512_xor_si512(D, _mm512_xor_si512(_mm512_and_si512(b, _mm512_srli_epi64(c, 2)),
                         _mm512_and_si512(_mm512_xor_si512(a, b), _mm512_srli_epi64(d, 2))));
    a = A; b = B; c = C; d = D;
    A = _mm512_xor_si512(_mm512_and_si512(a, _mm512_srli_epi64(a, 4)),
                         _mm512_and_si512(b, _mm512_srli_epi64(b, 4)));
    B = _mm512_xor_si512(_mm512_and_si512(a, _mm512_srli_epi64(b, 4)),
                         _mm512_and_si512(b, _mm512_srli_epi64(_mm512_xor_si512(a, b), 4)));
    C = _mm512_xor_si512(C, _mm512_xor_si512(_mm512_and_si512(a, _mm512_srli_epi64(c, 4)),
                         _mm512_and_si512(b, _mm512_srli_epi64(d, 4))));
    D = _mm512_xor_si512(D, _mm512_xor_si512(_mm512_and_si512(b, _mm512_srli_epi64(c, 4)),
                         _mm512_and_si512(_mm512_xor_si512(a, b), _mm512_srli_epi64(d, 4))));
    a = A; b = B; c = C; d = D;
    A = _mm512_xor_si512(_mm512_and_si512(a, _mm512_srli_epi64(a, 8)),
                         _mm512_and_si512(b, _mm512_srli_epi64(b, 8)));
    B = _mm512_xor_si512(_mm512_and_si512(a, _mm512_srli_epi64(b, 8)),
                         _mm512_and_si512(b, _mm512_srli_epi64(_mm512_xor_si512(a, b), 8)));
    C = _mm512_xor_si512(C, _mm512_xor_si512(_mm512_and_si512(a, _mm512_srli_epi64(c, 8)),
                         _mm512_and_si512(b, _mm512_srli_epi64(d, 8))));
    D = _mm512_xor_si512(D, _mm512_xor_si512(_mm512_and_si512(b, _mm512_srli_epi64(c, 8)),
                         _mm512_and_si512(_mm512_xor_si512(a, b), _mm512_srli_epi64(d, 8))));
    a = A; b = B; c = C; d = D;
    C = _mm512_xor_si512(C, _mm512_xor_si512(_mm512_and_si512(a, _mm512_srli_epi64(c, 16)),
                         _mm512_and_si512(b, _mm512_srli_epi64(d, 16))));
    D = _mm512_xor_si512(D, _mm512_xor_si512(_mm512_and_si512(b, _mm512_srli_epi64(c, 16)),
                         _mm512_and_si512(_mm512_xor_si512(a, b), _mm512_srli_epi64(d, 16))));
    a = _mm512_xor_si512(C, _mm512_srli_epi64(C, 1));
    b = _mm512_xor_si512(D, _mm512_srli_epi64(D, 1));
    i0 = _mm512_xor_si512(X, Y);
    i1 = _mm512_or_si512(b, _mm512_xor_si512(M, _mm512_or_si512(i0, a)));
    return _mm512_srlv_epi64(_mm512_or_si512(_mm512_slli_epi64(curve__part1by1_avx512_64(i1), 1),
                                             curve__part1by1_avx512_64(i0)),
                             _mm512_set1_epi64((long long)(64 - 2 * bits)));
}
CURVE_ATTR_AVX512 CURVE_INLINE void curve__hilbert2_64_d2xy_avx512(__m512i di, int bits,
                                                __m512i *CURVE_RESTRICT x,
                                                __m512i *CURVE_RESTRICT y)
{
    const __m512i M = _mm512_set1_epi64((long long)0xFFFFFFFFull);
    __m512i i, i0, i1, t0, t1, a;
    i = _mm512_sllv_epi64(di, _mm512_set1_epi64((long long)(64 - 2 * bits)));
    i0 = curve__compact1by1_avx512_64(i);
    i1 = curve__compact1by1_avx512_64(_mm512_srli_epi64(i, 1));
    t0 = _mm512_xor_si512(_mm512_or_si512(i0, i1), M);
    t1 = _mm512_and_si512(i0, i1);
    a = _mm512_or_si512(_mm512_and_si512(_mm512_xor_si512(i0, M), curve__prefixscan_avx512_64(t1)),
                        _mm512_and_si512(i0, curve__prefixscan_avx512_64(t0)));
    *x = _mm512_srlv_epi64(_mm512_xor_si512(a, i1),
                           _mm512_set1_epi64((long long)(32 - bits)));
    *y = _mm512_srlv_epi64(_mm512_xor_si512(_mm512_xor_si512(a, i0), i1),
                           _mm512_set1_epi64((long long)(32 - bits)));
}
#endif /* CURVE_HAVE_AVX512 */

/* Named complete-batch kernels. Public curve_*_batch functions pick the
 * compile-time widest ISA. curve_dispatch_init selects among compiled
 * kernels that the CPU actually runs. 64-bit Morton/Hilbert have no NEON
 * body (scalar on aarch64). Gray decode is scalar. */

CURVE_INLINE void curve__gray_encode_batch_scalar(
    const uint64_t *CURVE_RESTRICT x, uint64_t *CURVE_RESTRICT out, size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) out[i] = curve_gray_encode(x[i]);
}

#if defined(CURVE_HAVE_NEON)
CURVE_INLINE void curve__gray_encode_batch_neon(
    const uint64_t *CURVE_RESTRICT x, uint64_t *CURVE_RESTRICT out, size_t n)
{
    size_t i = 0;
    for (; i < n - (n & 1u); i += 2) {
        uint64x2_t v = vld1q_u64(x + i);
        vst1q_u64(out + i, veorq_u64(v, vshrq_n_u64(v, 1)));
    }
    for (; i < n; ++i) out[i] = curve_gray_encode(x[i]);
}
#endif

#if defined(CURVE_HAVE_AVX2) || defined(CURVE_DISPATCH_X86)
CURVE_ATTR_AVX2 CURVE_INLINE void curve__gray_encode_batch_avx2(
    const uint64_t *CURVE_RESTRICT x, uint64_t *CURVE_RESTRICT out, size_t n)
{
    size_t i = 0;
    for (; i < n - (n & 3u); i += 4) {
        __m256i v = _mm256_loadu_si256((const __m256i *)(const void *)(x + i));
        _mm256_storeu_si256((__m256i *)(void *)(out + i),
                            _mm256_xor_si256(v, _mm256_srli_epi64(v, 1)));
    }
    for (; i < n; ++i) out[i] = curve_gray_encode(x[i]);
}
#endif

#if defined(CURVE_HAVE_AVX512) || defined(CURVE_DISPATCH_X86)
CURVE_ATTR_AVX512 CURVE_INLINE void curve__gray_encode_batch_avx512(
    const uint64_t *CURVE_RESTRICT x, uint64_t *CURVE_RESTRICT out, size_t n)
{
    size_t i = 0;
    for (; i < n - (n & 7u); i += 8) {
        __m512i v = _mm512_loadu_si512((const void *)(x + i));
        _mm512_storeu_si512((void *)(out + i),
                            _mm512_xor_si512(v, _mm512_srli_epi64(v, 1)));
    }
    for (; i < n; ++i) out[i] = curve_gray_encode(x[i]);
}
#endif

CURVE_INLINE void curve__morton2_batch_scalar(
    const uint32_t *CURVE_RESTRICT x, const uint32_t *CURVE_RESTRICT y,
    uint32_t *CURVE_RESTRICT out, size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) out[i] = curve_morton2(x[i], y[i]);
}

CURVE_INLINE void curve__morton2_decode_batch_scalar(
    const uint32_t *CURVE_RESTRICT d, uint32_t *CURVE_RESTRICT x,
    uint32_t *CURVE_RESTRICT y, size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) curve_morton2_decode(d[i], x + i, y + i);
}

#if defined(CURVE_HAVE_NEON)
CURVE_INLINE void curve__morton2_batch_neon(
    const uint32_t *CURVE_RESTRICT x, const uint32_t *CURVE_RESTRICT y,
    uint32_t *CURVE_RESTRICT out, size_t n)
{
    size_t i = 0;
    for (; i < n - (n & 3u); i += 4) {
        uint32x4_t vx = vld1q_u32(x + i);
        uint32x4_t vy = vld1q_u32(y + i);
        vst1q_u32(out + i, vorrq_u32(curve__part1by1_u32x4(vx),
                                     vshlq_n_u32(curve__part1by1_u32x4(vy), 1)));
    }
    for (; i < n; ++i) out[i] = curve_morton2(x[i], y[i]);
}
CURVE_INLINE void curve__morton2_decode_batch_neon(
    const uint32_t *CURVE_RESTRICT d, uint32_t *CURVE_RESTRICT x,
    uint32_t *CURVE_RESTRICT y, size_t n)
{
    size_t i = 0;
    for (; i < n - (n & 3u); i += 4) {
        uint32x4_t vd = vld1q_u32(d + i);
        vst1q_u32(x + i, curve__compact1by1_u32x4(vd));
        vst1q_u32(y + i, curve__compact1by1_u32x4(vshrq_n_u32(vd, 1)));
    }
    for (; i < n; ++i) curve_morton2_decode(d[i], x + i, y + i);
}
#endif

#if defined(CURVE_HAVE_AVX2) || defined(CURVE_DISPATCH_X86)
CURVE_ATTR_AVX2 CURVE_INLINE void curve__morton2_batch_avx2(
    const uint32_t *CURVE_RESTRICT x, const uint32_t *CURVE_RESTRICT y,
    uint32_t *CURVE_RESTRICT out, size_t n)
{
    size_t i = 0;
    for (; i < n - (n & 7u); i += 8) {
        __m256i vx = _mm256_loadu_si256((const __m256i *)(const void *)(x + i));
        __m256i vy = _mm256_loadu_si256((const __m256i *)(const void *)(y + i));
        _mm256_storeu_si256((__m256i *)(void *)(out + i),
            _mm256_or_si256(curve__part1by1_avx2(vx),
                            _mm256_slli_epi32(curve__part1by1_avx2(vy), 1)));
    }
    for (; i < n; ++i) out[i] = curve_morton2(x[i], y[i]);
}
CURVE_ATTR_AVX2 CURVE_INLINE void curve__morton2_decode_batch_avx2(
    const uint32_t *CURVE_RESTRICT d, uint32_t *CURVE_RESTRICT x,
    uint32_t *CURVE_RESTRICT y, size_t n)
{
    size_t i = 0;
    for (; i < n - (n & 7u); i += 8) {
        __m256i vd = _mm256_loadu_si256((const __m256i *)(const void *)(d + i));
        _mm256_storeu_si256((__m256i *)(void *)(x + i), curve__compact1by1_avx2(vd));
        _mm256_storeu_si256((__m256i *)(void *)(y + i),
                            curve__compact1by1_avx2(_mm256_srli_epi32(vd, 1)));
    }
    for (; i < n; ++i) curve_morton2_decode(d[i], x + i, y + i);
}
#endif

#if defined(CURVE_HAVE_AVX512) || defined(CURVE_DISPATCH_X86)
CURVE_ATTR_AVX512 CURVE_INLINE void curve__morton2_batch_avx512(
    const uint32_t *CURVE_RESTRICT x, const uint32_t *CURVE_RESTRICT y,
    uint32_t *CURVE_RESTRICT out, size_t n)
{
    size_t i = 0;
    for (; i < n - (n & 15u); i += 16) {
        __m512i vx = _mm512_loadu_si512((const void *)(x + i));
        __m512i vy = _mm512_loadu_si512((const void *)(y + i));
        _mm512_storeu_si512((void *)(out + i),
            _mm512_or_si512(curve__part1by1_avx512(vx),
                            _mm512_slli_epi32(curve__part1by1_avx512(vy), 1)));
    }
    for (; i < n; ++i) out[i] = curve_morton2(x[i], y[i]);
}
CURVE_ATTR_AVX512 CURVE_INLINE void curve__morton2_decode_batch_avx512(
    const uint32_t *CURVE_RESTRICT d, uint32_t *CURVE_RESTRICT x,
    uint32_t *CURVE_RESTRICT y, size_t n)
{
    size_t i = 0;
    for (; i < n - (n & 15u); i += 16) {
        __m512i vd = _mm512_loadu_si512((const void *)(d + i));
        _mm512_storeu_si512((void *)(x + i), curve__compact1by1_avx512(vd));
        _mm512_storeu_si512((void *)(y + i),
                            curve__compact1by1_avx512(_mm512_srli_epi32(vd, 1)));
    }
    for (; i < n; ++i) curve_morton2_decode(d[i], x + i, y + i);
}
#endif

CURVE_INLINE void curve__morton2_64_batch_scalar(
    const uint32_t *CURVE_RESTRICT x, const uint32_t *CURVE_RESTRICT y,
    uint64_t *CURVE_RESTRICT out, size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) out[i] = curve_morton2_64(x[i], y[i]);
}

CURVE_INLINE void curve__morton2_64_decode_batch_scalar(
    const uint64_t *CURVE_RESTRICT d, uint32_t *CURVE_RESTRICT x,
    uint32_t *CURVE_RESTRICT y, size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) curve_morton2_64_decode(d[i], x + i, y + i);
}

CURVE_INLINE void curve__morton3_batch_scalar(
    const uint32_t *CURVE_RESTRICT x, const uint32_t *CURVE_RESTRICT y,
    const uint32_t *CURVE_RESTRICT z, uint64_t *CURVE_RESTRICT out, size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) out[i] = curve_morton3(x[i], y[i], z[i]);
}

CURVE_INLINE void curve__morton3_decode_batch_scalar(
    const uint64_t *CURVE_RESTRICT d, uint32_t *CURVE_RESTRICT x,
    uint32_t *CURVE_RESTRICT y, uint32_t *CURVE_RESTRICT z, size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) curve_morton3_decode(d[i], x + i, y + i, z + i);
}

#if defined(CURVE_HAVE_AVX2) || defined(CURVE_DISPATCH_X86)
CURVE_ATTR_AVX2 CURVE_INLINE void curve__morton2_64_batch_avx2(
    const uint32_t *CURVE_RESTRICT x, const uint32_t *CURVE_RESTRICT y,
    uint64_t *CURVE_RESTRICT out, size_t n)
{
    size_t i = 0;
    for (; i < n - (n & 3u); i += 4) {
        __m128i x32 = _mm_loadu_si128((const __m128i *)(const void *)(x + i));
        __m128i y32 = _mm_loadu_si128((const __m128i *)(const void *)(y + i));
        __m256i vx = curve__part1by1_avx2_64(_mm256_cvtepu32_epi64(x32));
        __m256i vy = curve__part1by1_avx2_64(_mm256_cvtepu32_epi64(y32));
        _mm256_storeu_si256((__m256i *)(void *)(out + i),
                            _mm256_or_si256(vx, _mm256_slli_epi64(vy, 1)));
    }
    for (; i < n; ++i) out[i] = curve_morton2_64(x[i], y[i]);
}
CURVE_ATTR_AVX2 CURVE_INLINE void curve__morton2_64_decode_batch_avx2(
    const uint64_t *CURVE_RESTRICT d, uint32_t *CURVE_RESTRICT x,
    uint32_t *CURVE_RESTRICT y, size_t n)
{
    size_t i = 0;
    for (; i < n - (n & 3u); i += 4) {
        __m256i vd = _mm256_loadu_si256((const __m256i *)(const void *)(d + i));
        _mm_storeu_si128((__m128i *)(void *)(x + i),
            curve__pack64lo32_avx2(curve__compact1by1_avx2_64(vd)));
        _mm_storeu_si128((__m128i *)(void *)(y + i),
            curve__pack64lo32_avx2(curve__compact1by1_avx2_64(_mm256_srli_epi64(vd, 1))));
    }
    for (; i < n; ++i) curve_morton2_64_decode(d[i], x + i, y + i);
}
CURVE_ATTR_AVX2 CURVE_INLINE void curve__morton3_batch_avx2(
    const uint32_t *CURVE_RESTRICT x, const uint32_t *CURVE_RESTRICT y,
    const uint32_t *CURVE_RESTRICT z, uint64_t *CURVE_RESTRICT out, size_t n)
{
    size_t i = 0;
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
    for (; i < n; ++i) out[i] = curve_morton3(x[i], y[i], z[i]);
}
CURVE_ATTR_AVX2 CURVE_INLINE void curve__morton3_decode_batch_avx2(
    const uint64_t *CURVE_RESTRICT d, uint32_t *CURVE_RESTRICT x,
    uint32_t *CURVE_RESTRICT y, uint32_t *CURVE_RESTRICT z, size_t n)
{
    size_t i = 0;
    for (; i < n - (n & 3u); i += 4) {
        __m256i vd = _mm256_loadu_si256((const __m256i *)(const void *)(d + i));
        _mm_storeu_si128((__m128i *)(void *)(x + i),
            curve__pack64lo32_avx2(curve__compact1by2_avx2_64(vd)));
        _mm_storeu_si128((__m128i *)(void *)(y + i),
            curve__pack64lo32_avx2(curve__compact1by2_avx2_64(_mm256_srli_epi64(vd, 1))));
        _mm_storeu_si128((__m128i *)(void *)(z + i),
            curve__pack64lo32_avx2(curve__compact1by2_avx2_64(_mm256_srli_epi64(vd, 2))));
    }
    for (; i < n; ++i) curve_morton3_decode(d[i], x + i, y + i, z + i);
}
#endif

#if defined(CURVE_HAVE_AVX512) || defined(CURVE_DISPATCH_X86)
CURVE_ATTR_AVX512 CURVE_INLINE void curve__morton2_64_batch_avx512(
    const uint32_t *CURVE_RESTRICT x, const uint32_t *CURVE_RESTRICT y,
    uint64_t *CURVE_RESTRICT out, size_t n)
{
    size_t i = 0;
    for (; i < n - (n & 7u); i += 8) {
        __m256i x32 = _mm256_loadu_si256((const __m256i *)(const void *)(x + i));
        __m256i y32 = _mm256_loadu_si256((const __m256i *)(const void *)(y + i));
        __m512i vx = curve__part1by1_avx512_64(_mm512_cvtepu32_epi64(x32));
        __m512i vy = curve__part1by1_avx512_64(_mm512_cvtepu32_epi64(y32));
        _mm512_storeu_si512((void *)(out + i),
                            _mm512_or_si512(vx, _mm512_slli_epi64(vy, 1)));
    }
    for (; i < n; ++i) out[i] = curve_morton2_64(x[i], y[i]);
}
CURVE_ATTR_AVX512 CURVE_INLINE void curve__morton2_64_decode_batch_avx512(
    const uint64_t *CURVE_RESTRICT d, uint32_t *CURVE_RESTRICT x,
    uint32_t *CURVE_RESTRICT y, size_t n)
{
    size_t i = 0;
    for (; i < n - (n & 7u); i += 8) {
        __m512i vd = _mm512_loadu_si512((const void *)(d + i));
        _mm256_storeu_si256((__m256i *)(void *)(x + i),
            _mm512_cvtepi64_epi32(curve__compact1by1_avx512_64(vd)));
        _mm256_storeu_si256((__m256i *)(void *)(y + i),
            _mm512_cvtepi64_epi32(curve__compact1by1_avx512_64(_mm512_srli_epi64(vd, 1))));
    }
    for (; i < n; ++i) curve_morton2_64_decode(d[i], x + i, y + i);
}
CURVE_ATTR_AVX512 CURVE_INLINE void curve__morton3_batch_avx512(
    const uint32_t *CURVE_RESTRICT x, const uint32_t *CURVE_RESTRICT y,
    const uint32_t *CURVE_RESTRICT z, uint64_t *CURVE_RESTRICT out, size_t n)
{
    size_t i = 0;
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
    for (; i < n; ++i) out[i] = curve_morton3(x[i], y[i], z[i]);
}
CURVE_ATTR_AVX512 CURVE_INLINE void curve__morton3_decode_batch_avx512(
    const uint64_t *CURVE_RESTRICT d, uint32_t *CURVE_RESTRICT x,
    uint32_t *CURVE_RESTRICT y, uint32_t *CURVE_RESTRICT z, size_t n)
{
    size_t i = 0;
    for (; i < n - (n & 7u); i += 8) {
        __m512i vd = _mm512_loadu_si512((const void *)(d + i));
        _mm256_storeu_si256((__m256i *)(void *)(x + i),
            _mm512_cvtepi64_epi32(curve__compact1by2_avx512_64(vd)));
        _mm256_storeu_si256((__m256i *)(void *)(y + i),
            _mm512_cvtepi64_epi32(curve__compact1by2_avx512_64(_mm512_srli_epi64(vd, 1))));
        _mm256_storeu_si256((__m256i *)(void *)(z + i),
            _mm512_cvtepi64_epi32(curve__compact1by2_avx512_64(_mm512_srli_epi64(vd, 2))));
    }
    for (; i < n; ++i) curve_morton3_decode(d[i], x + i, y + i, z + i);
}
#endif

CURVE_INLINE void curve__hilbert_xy2d_batch_scalar(
    const uint32_t *CURVE_RESTRICT x, const uint32_t *CURVE_RESTRICT y,
    int bits, uint32_t *CURVE_RESTRICT out, size_t n)
{
    size_t i;
    bits = curve__clamp_bits2(bits);
    for (i = 0; i < n; ++i) out[i] = curve_hilbert_xy2d(x[i], y[i], bits);
}

CURVE_INLINE void curve__hilbert_d2xy_batch_scalar(
    const uint32_t *CURVE_RESTRICT d, int bits,
    uint32_t *CURVE_RESTRICT x, uint32_t *CURVE_RESTRICT y, size_t n)
{
    size_t i;
    bits = curve__clamp_bits2(bits);
    for (i = 0; i < n; ++i) curve_hilbert_d2xy(d[i], bits, x + i, y + i);
}

#if defined(CURVE_HAVE_NEON)
CURVE_INLINE void curve__hilbert_xy2d_batch_neon(
    const uint32_t *CURVE_RESTRICT x, const uint32_t *CURVE_RESTRICT y,
    int bits, uint32_t *CURVE_RESTRICT out, size_t n)
{
    size_t i = 0;
    bits = curve__clamp_bits2(bits);
    for (; i < n - (n & 3u); i += 4) {
        uint32x4_t vx = vld1q_u32(x + i);
        uint32x4_t vy = vld1q_u32(y + i);
        vst1q_u32(out + i, curve__hilbert_xy2d_u32x4(vx, vy, bits));
    }
    for (; i < n; ++i) out[i] = curve_hilbert_xy2d(x[i], y[i], bits);
}
CURVE_INLINE void curve__hilbert_d2xy_batch_neon(
    const uint32_t *CURVE_RESTRICT d, int bits,
    uint32_t *CURVE_RESTRICT x, uint32_t *CURVE_RESTRICT y, size_t n)
{
    size_t i = 0;
    bits = curve__clamp_bits2(bits);
    for (; i < n - (n & 3u); i += 4) {
        uint32x4_t vd = vld1q_u32(d + i);
        uint32x4_t vx, vy;
        curve__hilbert_d2xy_u32x4(vd, bits, &vx, &vy);
        vst1q_u32(x + i, vx);
        vst1q_u32(y + i, vy);
    }
    for (; i < n; ++i) curve_hilbert_d2xy(d[i], bits, x + i, y + i);
}
#endif

#if defined(CURVE_HAVE_AVX2) || defined(CURVE_DISPATCH_X86)
CURVE_ATTR_AVX2 CURVE_INLINE void curve__hilbert_xy2d_batch_avx2(
    const uint32_t *CURVE_RESTRICT x, const uint32_t *CURVE_RESTRICT y,
    int bits, uint32_t *CURVE_RESTRICT out, size_t n)
{
    size_t i = 0;
    bits = curve__clamp_bits2(bits);
    for (; i < n - (n & 7u); i += 8) {
        __m256i vx = _mm256_loadu_si256((const __m256i *)(const void *)(x + i));
        __m256i vy = _mm256_loadu_si256((const __m256i *)(const void *)(y + i));
        _mm256_storeu_si256((__m256i *)(void *)(out + i),
                            curve__hilbert_xy2d_avx2(vx, vy, bits));
    }
    for (; i < n; ++i) out[i] = curve_hilbert_xy2d(x[i], y[i], bits);
}
CURVE_ATTR_AVX2 CURVE_INLINE void curve__hilbert_d2xy_batch_avx2(
    const uint32_t *CURVE_RESTRICT d, int bits,
    uint32_t *CURVE_RESTRICT x, uint32_t *CURVE_RESTRICT y, size_t n)
{
    size_t i = 0;
    bits = curve__clamp_bits2(bits);
    for (; i < n - (n & 7u); i += 8) {
        __m256i vd = _mm256_loadu_si256((const __m256i *)(const void *)(d + i));
        __m256i vx, vy;
        curve__hilbert_d2xy_avx2(vd, bits, &vx, &vy);
        _mm256_storeu_si256((__m256i *)(void *)(x + i), vx);
        _mm256_storeu_si256((__m256i *)(void *)(y + i), vy);
    }
    for (; i < n; ++i) curve_hilbert_d2xy(d[i], bits, x + i, y + i);
}
#endif

#if defined(CURVE_HAVE_AVX512) || defined(CURVE_DISPATCH_X86)
CURVE_ATTR_AVX512 CURVE_INLINE void curve__hilbert_xy2d_batch_avx512(
    const uint32_t *CURVE_RESTRICT x, const uint32_t *CURVE_RESTRICT y,
    int bits, uint32_t *CURVE_RESTRICT out, size_t n)
{
    size_t i = 0;
    bits = curve__clamp_bits2(bits);
    for (; i < n - (n & 15u); i += 16) {
        __m512i vx = _mm512_loadu_si512((const void *)(x + i));
        __m512i vy = _mm512_loadu_si512((const void *)(y + i));
        _mm512_storeu_si512((void *)(out + i),
                            curve__hilbert_xy2d_avx512(vx, vy, bits));
    }
    for (; i < n; ++i) out[i] = curve_hilbert_xy2d(x[i], y[i], bits);
}
CURVE_ATTR_AVX512 CURVE_INLINE void curve__hilbert_d2xy_batch_avx512(
    const uint32_t *CURVE_RESTRICT d, int bits,
    uint32_t *CURVE_RESTRICT x, uint32_t *CURVE_RESTRICT y, size_t n)
{
    size_t i = 0;
    bits = curve__clamp_bits2(bits);
    for (; i < n - (n & 15u); i += 16) {
        __m512i vd = _mm512_loadu_si512((const void *)(d + i));
        __m512i vx, vy;
        curve__hilbert_d2xy_avx512(vd, bits, &vx, &vy);
        _mm512_storeu_si512((void *)(x + i), vx);
        _mm512_storeu_si512((void *)(y + i), vy);
    }
    for (; i < n; ++i) curve_hilbert_d2xy(d[i], bits, x + i, y + i);
}
#endif

CURVE_INLINE void curve__hilbert2_64_xy2d_batch_scalar(
    const uint32_t *CURVE_RESTRICT x, const uint32_t *CURVE_RESTRICT y,
    int bits, uint64_t *CURVE_RESTRICT out, size_t n)
{
    size_t i;
    bits = curve__clamp_bits2_64(bits);
    for (i = 0; i < n; ++i) out[i] = curve_hilbert2_64_xy2d(x[i], y[i], bits);
}

CURVE_INLINE void curve__hilbert2_64_d2xy_batch_scalar(
    const uint64_t *CURVE_RESTRICT d, int bits,
    uint32_t *CURVE_RESTRICT x, uint32_t *CURVE_RESTRICT y, size_t n)
{
    size_t i;
    bits = curve__clamp_bits2_64(bits);
    for (i = 0; i < n; ++i) curve_hilbert2_64_d2xy(d[i], bits, x + i, y + i);
}

#if defined(CURVE_HAVE_AVX2) || defined(CURVE_DISPATCH_X86)
CURVE_ATTR_AVX2 CURVE_INLINE void curve__hilbert2_64_xy2d_batch_avx2(
    const uint32_t *CURVE_RESTRICT x, const uint32_t *CURVE_RESTRICT y,
    int bits, uint64_t *CURVE_RESTRICT out, size_t n)
{
    size_t i = 0;
    bits = curve__clamp_bits2_64(bits);
    for (; i < n - (n & 3u); i += 4) {
        __m128i x32 = _mm_loadu_si128((const __m128i *)(const void *)(x + i));
        __m128i y32 = _mm_loadu_si128((const __m128i *)(const void *)(y + i));
        __m256i vx = _mm256_cvtepu32_epi64(x32);
        __m256i vy = _mm256_cvtepu32_epi64(y32);
        _mm256_storeu_si256((__m256i *)(void *)(out + i),
                            curve__hilbert2_64_xy2d_avx2(vx, vy, bits));
    }
    for (; i < n; ++i) out[i] = curve_hilbert2_64_xy2d(x[i], y[i], bits);
}
CURVE_ATTR_AVX2 CURVE_INLINE void curve__hilbert2_64_d2xy_batch_avx2(
    const uint64_t *CURVE_RESTRICT d, int bits,
    uint32_t *CURVE_RESTRICT x, uint32_t *CURVE_RESTRICT y, size_t n)
{
    size_t i = 0;
    bits = curve__clamp_bits2_64(bits);
    for (; i < n - (n & 3u); i += 4) {
        __m256i vd = _mm256_loadu_si256((const __m256i *)(const void *)(d + i));
        __m256i vx, vy;
        curve__hilbert2_64_d2xy_avx2(vd, bits, &vx, &vy);
        _mm_storeu_si128((__m128i *)(void *)(x + i), curve__pack64lo32_avx2(vx));
        _mm_storeu_si128((__m128i *)(void *)(y + i), curve__pack64lo32_avx2(vy));
    }
    for (; i < n; ++i) curve_hilbert2_64_d2xy(d[i], bits, x + i, y + i);
}
#endif

#if defined(CURVE_HAVE_AVX512) || defined(CURVE_DISPATCH_X86)
CURVE_ATTR_AVX512 CURVE_INLINE void curve__hilbert2_64_xy2d_batch_avx512(
    const uint32_t *CURVE_RESTRICT x, const uint32_t *CURVE_RESTRICT y,
    int bits, uint64_t *CURVE_RESTRICT out, size_t n)
{
    size_t i = 0;
    bits = curve__clamp_bits2_64(bits);
    for (; i < n - (n & 7u); i += 8) {
        __m256i x32 = _mm256_loadu_si256((const __m256i *)(const void *)(x + i));
        __m256i y32 = _mm256_loadu_si256((const __m256i *)(const void *)(y + i));
        __m512i vx = _mm512_cvtepu32_epi64(x32);
        __m512i vy = _mm512_cvtepu32_epi64(y32);
        _mm512_storeu_si512((void *)(out + i),
                            curve__hilbert2_64_xy2d_avx512(vx, vy, bits));
    }
    for (; i < n; ++i) out[i] = curve_hilbert2_64_xy2d(x[i], y[i], bits);
}
CURVE_ATTR_AVX512 CURVE_INLINE void curve__hilbert2_64_d2xy_batch_avx512(
    const uint64_t *CURVE_RESTRICT d, int bits,
    uint32_t *CURVE_RESTRICT x, uint32_t *CURVE_RESTRICT y, size_t n)
{
    size_t i = 0;
    bits = curve__clamp_bits2_64(bits);
    for (; i < n - (n & 7u); i += 8) {
        __m512i vd = _mm512_loadu_si512((const void *)(d + i));
        __m512i vx, vy;
        curve__hilbert2_64_d2xy_avx512(vd, bits, &vx, &vy);
        _mm256_storeu_si256((__m256i *)(void *)(x + i), _mm512_cvtepi64_epi32(vx));
        _mm256_storeu_si256((__m256i *)(void *)(y + i), _mm512_cvtepi64_epi32(vy));
    }
    for (; i < n; ++i) curve_hilbert2_64_d2xy(d[i], bits, x + i, y + i);
}
#endif

#if defined(CURVE_HAVE_AVX512)
#  define CURVE__BATCH_SIMD(fn) curve__##fn##_avx512
#  define CURVE__BATCH_AVX(fn)  curve__##fn##_avx512
#elif defined(CURVE_HAVE_AVX2)
#  define CURVE__BATCH_SIMD(fn) curve__##fn##_avx2
#  define CURVE__BATCH_AVX(fn)  curve__##fn##_avx2
#elif defined(CURVE_HAVE_NEON)
#  define CURVE__BATCH_SIMD(fn) curve__##fn##_neon
#  define CURVE__BATCH_AVX(fn)  curve__##fn##_scalar
#else
#  define CURVE__BATCH_SIMD(fn) curve__##fn##_scalar
#  define CURVE__BATCH_AVX(fn)  curve__##fn##_scalar
#endif

CURVE_INLINE void curve_gray_encode_batch(const uint64_t *CURVE_RESTRICT x,
                                          uint64_t *CURVE_RESTRICT out, size_t n)
{
    CURVE__BATCH_SIMD(gray_encode_batch)(x, out, n);
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
    CURVE__BATCH_SIMD(morton2_batch)(x, y, out, n);
}

CURVE_INLINE void curve_morton2_decode_batch(const uint32_t *CURVE_RESTRICT d,
                                             uint32_t *CURVE_RESTRICT x,
                                             uint32_t *CURVE_RESTRICT y,
                                             size_t n)
{
    CURVE__BATCH_SIMD(morton2_decode_batch)(d, x, y, n);
}

CURVE_INLINE void curve_morton2_64_batch(const uint32_t *CURVE_RESTRICT x,
                                         const uint32_t *CURVE_RESTRICT y,
                                         uint64_t *CURVE_RESTRICT out, size_t n)
{
    CURVE__BATCH_AVX(morton2_64_batch)(x, y, out, n);
}

CURVE_INLINE void curve_morton2_64_decode_batch(const uint64_t *CURVE_RESTRICT d,
                                                uint32_t *CURVE_RESTRICT x,
                                                uint32_t *CURVE_RESTRICT y,
                                                size_t n)
{
    CURVE__BATCH_AVX(morton2_64_decode_batch)(d, x, y, n);
}

CURVE_INLINE void curve_morton3_batch(const uint32_t *CURVE_RESTRICT x,
                                      const uint32_t *CURVE_RESTRICT y,
                                      const uint32_t *CURVE_RESTRICT z,
                                      uint64_t *CURVE_RESTRICT out, size_t n)
{
    CURVE__BATCH_AVX(morton3_batch)(x, y, z, out, n);
}

CURVE_INLINE void curve_morton3_decode_batch(const uint64_t *CURVE_RESTRICT d,
                                             uint32_t *CURVE_RESTRICT x,
                                             uint32_t *CURVE_RESTRICT y,
                                             uint32_t *CURVE_RESTRICT z,
                                             size_t n)
{
    CURVE__BATCH_AVX(morton3_decode_batch)(d, x, y, z, n);
}

CURVE_INLINE void curve_hilbert_xy2d_batch(const uint32_t *CURVE_RESTRICT x,
                                           const uint32_t *CURVE_RESTRICT y,
                                           int bits, uint32_t *CURVE_RESTRICT out,
                                           size_t n)
{
    CURVE__BATCH_SIMD(hilbert_xy2d_batch)(x, y, bits, out, n);
}

CURVE_INLINE void curve_hilbert_d2xy_batch(const uint32_t *CURVE_RESTRICT d,
                                           int bits,
                                           uint32_t *CURVE_RESTRICT x,
                                           uint32_t *CURVE_RESTRICT y,
                                           size_t n)
{
    CURVE__BATCH_SIMD(hilbert_d2xy_batch)(d, bits, x, y, n);
}

CURVE_INLINE void curve_hilbert2_64_xy2d_batch(const uint32_t *CURVE_RESTRICT x,
                                               const uint32_t *CURVE_RESTRICT y,
                                               int bits,
                                               uint64_t *CURVE_RESTRICT out,
                                               size_t n)
{
    CURVE__BATCH_AVX(hilbert2_64_xy2d_batch)(x, y, bits, out, n);
}

CURVE_INLINE void curve_hilbert2_64_d2xy_batch(const uint64_t *CURVE_RESTRICT d,
                                               int bits,
                                               uint32_t *CURVE_RESTRICT x,
                                               uint32_t *CURVE_RESTRICT y,
                                               size_t n)
{
    CURVE__BATCH_AVX(hilbert2_64_d2xy_batch)(d, bits, x, y, n);
}

#undef CURVE__BATCH_SIMD
#undef CURVE__BATCH_AVX

/* -------------------------------------------------------------------------- */
/* Runtime ISA dispatch (caller-owned; no process globals)                    */
/*                                                                            */
/* Default public batch APIs stay compile-time.  curve_dispatch_init fills a  */
/* vtable from compiled ∩ runtime ∩ allow.  allow == 0 means "any compiled    */
/* kernel the CPU can run".  SCALAR is always available.                      */
/*                                                                            */
/* CURVE_RUNTIME_DISPATCH (x86 gcc/clang) compiles extra AVX2/AVX-512         */
/* kernels with target attributes so a baseline binary can still pick them    */
/* at run time.  MSVC has no target attributes; the vtable only lists the     */
/* compile-time path.  aarch64 NEON is baseline, so no getauxval.             */
/* -------------------------------------------------------------------------- */

CURVE_INLINE unsigned curve_isa_compiled_features(void)
{
    unsigned f = CURVE_ISA_SCALAR;
#if defined(CURVE_HAVE_NEON)
    f |= CURVE_ISA_NEON;
#endif
#if defined(CURVE_HAVE_AVX2) || defined(CURVE_DISPATCH_X86)
    f |= CURVE_ISA_AVX2;
#endif
#if defined(CURVE_HAVE_AVX512) || defined(CURVE_DISPATCH_X86)
    f |= CURVE_ISA_AVX512;
#endif
#if defined(CURVE_HAVE_BMI2)
    f |= CURVE_ISA_BMI2;
#endif
    return f;
}

#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
CURVE_INLINE unsigned curve__x86_runtime_features(void)
{
    unsigned a1 = 0, b1 = 0, c1 = 0, d1 = 0;
    unsigned a7 = 0, b7 = 0, c7 = 0, d7 = 0;
    unsigned xcr0, f = 0;
#if defined(_MSC_VER)
    int r1[4], r7[4];
    __cpuidex(r1, 1, 0);
    a1 = (unsigned)r1[0]; b1 = (unsigned)r1[1];
    c1 = (unsigned)r1[2]; d1 = (unsigned)r1[3];
    __cpuidex(r7, 7, 0);
    a7 = (unsigned)r7[0]; b7 = (unsigned)r7[1];
    c7 = (unsigned)r7[2]; d7 = (unsigned)r7[3];
#elif defined(__GNUC__) || defined(__clang__)
    if (!__get_cpuid(1, &a1, &b1, &c1, &d1))
        return 0;
    if (__get_cpuid_max(0, NULL) >= 7u)
        (void)__get_cpuid_count(7, 0, &a7, &b7, &c7, &d7);
#else
    return 0;
#endif
    (void)a1; (void)b1; (void)d1; (void)a7; (void)c7; (void)d7;
    if (b7 & (1u << 8)) f |= CURVE_ISA_BMI2;
    if ((c1 & ((1u << 27) | (1u << 28))) != ((1u << 27) | (1u << 28)))
        return f;
#if defined(_MSC_VER)
    xcr0 = (unsigned)_xgetbv(0);
#elif defined(__GNUC__) || defined(__clang__)
    {
        unsigned eax, edx;
#  if defined(__GNUC__)
#    pragma GCC diagnostic push
#    pragma GCC diagnostic ignored "-Wpedantic"
#  endif
        __asm__ __volatile__(".byte 0x0f, 0x01, 0xd0"
                             : "=a"(eax), "=d"(edx) : "c"(0));
#  if defined(__GNUC__)
#    pragma GCC diagnostic pop
#  endif
        xcr0 = eax;
    }
#else
    return f;
#endif
    if ((xcr0 & 0x6u) != 0x6u)
        return f;
    if (b7 & (1u << 5)) f |= CURVE_ISA_AVX2;
    if ((b7 & (1u << 16)) && ((xcr0 & 0xE6u) == 0xE6u))
        f |= CURVE_ISA_AVX512;
    return f;
}
#endif

CURVE_INLINE unsigned curve_isa_runtime_features(void)
{
    unsigned f = CURVE_ISA_SCALAR;
#if defined(__aarch64__) || defined(_M_ARM64)
    f |= CURVE_ISA_NEON;
#elif defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
    f |= curve__x86_runtime_features();
#elif defined(__ARM_NEON) || defined(__ARM_NEON__)
    f |= CURVE_ISA_NEON;
#endif
    return f;
}

CURVE_INLINE unsigned curve_isa_detect(unsigned allow)
{
    unsigned f = curve_isa_compiled_features() & curve_isa_runtime_features();
    if (allow) f &= allow;
    f |= CURVE_ISA_SCALAR;
    return f;
}

CURVE_INLINE unsigned curve__isa_pick(unsigned features)
{
    if (features & CURVE_ISA_AVX512) return CURVE_ISA_AVX512;
    if (features & CURVE_ISA_AVX2)   return CURVE_ISA_AVX2;
    if (features & CURVE_ISA_NEON)   return CURVE_ISA_NEON;
    return CURVE_ISA_SCALAR;
}

CURVE_INLINE const char *curve_isa_name(unsigned features)
{
    if (features & CURVE_ISA_AVX512) return "avx512";
    if (features & CURVE_ISA_AVX2)   return "avx2";
    if (features & CURVE_ISA_NEON)   return "neon";
    return "scalar";
}

typedef struct curve_dispatch {
    unsigned features;
    void (*gray_encode_batch)(const uint64_t *CURVE_RESTRICT,
                              uint64_t *CURVE_RESTRICT, size_t);
    void (*gray_decode_batch)(const uint64_t *CURVE_RESTRICT,
                              uint64_t *CURVE_RESTRICT, size_t);
    void (*morton2_batch)(const uint32_t *CURVE_RESTRICT,
                          const uint32_t *CURVE_RESTRICT,
                          uint32_t *CURVE_RESTRICT, size_t);
    void (*morton2_decode_batch)(const uint32_t *CURVE_RESTRICT,
                                 uint32_t *CURVE_RESTRICT,
                                 uint32_t *CURVE_RESTRICT, size_t);
    void (*morton2_64_batch)(const uint32_t *CURVE_RESTRICT,
                             const uint32_t *CURVE_RESTRICT,
                             uint64_t *CURVE_RESTRICT, size_t);
    void (*morton2_64_decode_batch)(const uint64_t *CURVE_RESTRICT,
                                    uint32_t *CURVE_RESTRICT,
                                    uint32_t *CURVE_RESTRICT, size_t);
    void (*morton3_batch)(const uint32_t *CURVE_RESTRICT,
                          const uint32_t *CURVE_RESTRICT,
                          const uint32_t *CURVE_RESTRICT,
                          uint64_t *CURVE_RESTRICT, size_t);
    void (*morton3_decode_batch)(const uint64_t *CURVE_RESTRICT,
                                 uint32_t *CURVE_RESTRICT,
                                 uint32_t *CURVE_RESTRICT,
                                 uint32_t *CURVE_RESTRICT, size_t);
    void (*hilbert_xy2d_batch)(const uint32_t *CURVE_RESTRICT,
                               const uint32_t *CURVE_RESTRICT, int,
                               uint32_t *CURVE_RESTRICT, size_t);
    void (*hilbert_d2xy_batch)(const uint32_t *CURVE_RESTRICT, int,
                               uint32_t *CURVE_RESTRICT,
                               uint32_t *CURVE_RESTRICT, size_t);
    void (*hilbert2_64_xy2d_batch)(const uint32_t *CURVE_RESTRICT,
                                   const uint32_t *CURVE_RESTRICT, int,
                                   uint64_t *CURVE_RESTRICT, size_t);
    void (*hilbert2_64_d2xy_batch)(const uint64_t *CURVE_RESTRICT, int,
                                   uint32_t *CURVE_RESTRICT,
                                   uint32_t *CURVE_RESTRICT, size_t);
} curve_dispatch;

CURVE_INLINE void curve__dispatch_fill_scalar(curve_dispatch *d)
{
    d->gray_encode_batch = curve__gray_encode_batch_scalar;
    d->gray_decode_batch = curve_gray_decode_batch;
    d->morton2_batch = curve__morton2_batch_scalar;
    d->morton2_decode_batch = curve__morton2_decode_batch_scalar;
    d->morton2_64_batch = curve__morton2_64_batch_scalar;
    d->morton2_64_decode_batch = curve__morton2_64_decode_batch_scalar;
    d->morton3_batch = curve__morton3_batch_scalar;
    d->morton3_decode_batch = curve__morton3_decode_batch_scalar;
    d->hilbert_xy2d_batch = curve__hilbert_xy2d_batch_scalar;
    d->hilbert_d2xy_batch = curve__hilbert_d2xy_batch_scalar;
    d->hilbert2_64_xy2d_batch = curve__hilbert2_64_xy2d_batch_scalar;
    d->hilbert2_64_d2xy_batch = curve__hilbert2_64_d2xy_batch_scalar;
}

#if defined(CURVE_HAVE_NEON)
CURVE_INLINE void curve__dispatch_fill_neon(curve_dispatch *d)
{
    d->gray_encode_batch = curve__gray_encode_batch_neon;
    d->morton2_batch = curve__morton2_batch_neon;
    d->morton2_decode_batch = curve__morton2_decode_batch_neon;
    d->hilbert_xy2d_batch = curve__hilbert_xy2d_batch_neon;
    d->hilbert_d2xy_batch = curve__hilbert_d2xy_batch_neon;
}
#endif

#if defined(CURVE_HAVE_AVX2) || defined(CURVE_DISPATCH_X86)
CURVE_INLINE void curve__dispatch_fill_avx2(curve_dispatch *d)
{
    d->gray_encode_batch = curve__gray_encode_batch_avx2;
    d->morton2_batch = curve__morton2_batch_avx2;
    d->morton2_decode_batch = curve__morton2_decode_batch_avx2;
    d->morton2_64_batch = curve__morton2_64_batch_avx2;
    d->morton2_64_decode_batch = curve__morton2_64_decode_batch_avx2;
    d->morton3_batch = curve__morton3_batch_avx2;
    d->morton3_decode_batch = curve__morton3_decode_batch_avx2;
    d->hilbert_xy2d_batch = curve__hilbert_xy2d_batch_avx2;
    d->hilbert_d2xy_batch = curve__hilbert_d2xy_batch_avx2;
    d->hilbert2_64_xy2d_batch = curve__hilbert2_64_xy2d_batch_avx2;
    d->hilbert2_64_d2xy_batch = curve__hilbert2_64_d2xy_batch_avx2;
}
#endif

#if defined(CURVE_HAVE_AVX512) || defined(CURVE_DISPATCH_X86)
CURVE_INLINE void curve__dispatch_fill_avx512(curve_dispatch *d)
{
    d->gray_encode_batch = curve__gray_encode_batch_avx512;
    d->morton2_batch = curve__morton2_batch_avx512;
    d->morton2_decode_batch = curve__morton2_decode_batch_avx512;
    d->morton2_64_batch = curve__morton2_64_batch_avx512;
    d->morton2_64_decode_batch = curve__morton2_64_decode_batch_avx512;
    d->morton3_batch = curve__morton3_batch_avx512;
    d->morton3_decode_batch = curve__morton3_decode_batch_avx512;
    d->hilbert_xy2d_batch = curve__hilbert_xy2d_batch_avx512;
    d->hilbert_d2xy_batch = curve__hilbert_d2xy_batch_avx512;
    d->hilbert2_64_xy2d_batch = curve__hilbert2_64_xy2d_batch_avx512;
    d->hilbert2_64_d2xy_batch = curve__hilbert2_64_d2xy_batch_avx512;
}
#endif

CURVE_INLINE void curve_dispatch_init(curve_dispatch *d, unsigned allow)
{
    unsigned pick;
    if (!d) return;
    pick = curve__isa_pick(curve_isa_detect(allow));
    d->features = pick | CURVE_ISA_SCALAR;
    curve__dispatch_fill_scalar(d);
#if defined(CURVE_HAVE_AVX512) || defined(CURVE_DISPATCH_X86)
    if (pick == CURVE_ISA_AVX512) {
        curve__dispatch_fill_avx512(d);
        return;
    }
#endif
#if defined(CURVE_HAVE_AVX2) || defined(CURVE_DISPATCH_X86)
    if (pick == CURVE_ISA_AVX2) {
        curve__dispatch_fill_avx2(d);
        return;
    }
#endif
#if defined(CURVE_HAVE_NEON)
    if (pick == CURVE_ISA_NEON) {
        curve__dispatch_fill_neon(d);
        return;
    }
#endif
}

CURVE_INLINE const char *curve_selected_isa(const curve_dispatch *d)
{
    if (!d) return "none";
    return curve_isa_name(d->features);
}

#ifdef __cplusplus
}
#endif

#endif /* CURVES_H */
