/* Independent reference encoders for tests. Not part of the public API.
 *
 * Morton 2D/3D match the legacy even/odd-from-LSB layout of curve_morton2 /
 * curve_morton3. Wide n-D packing uses MSB-level-first, axis 0 first within
 * a level (the n-D Hilbert convention). */
#ifndef REF_CURVES_H
#define REF_CURVES_H

#include <stdint.h>
#include <stddef.h>

#if defined(__GNUC__) || defined(__clang__)
#  define REF_UNUSED __attribute__((unused))
#else
#  define REF_UNUSED
#endif

static REF_UNUSED uint32_t ref_morton2(uint32_t x, uint32_t y)
{
    uint32_t d = 0;
    int i;
    x &= 0xFFFFu;
    y &= 0xFFFFu;
    for (i = 0; i < 16; i++) {
        d |= ((x >> i) & 1u) << (2 * i);
        d |= ((y >> i) & 1u) << (2 * i + 1);
    }
    return d;
}

static REF_UNUSED void ref_morton2_decode(uint32_t d, uint32_t *x, uint32_t *y)
{
    uint32_t rx = 0, ry = 0;
    int i;
    for (i = 0; i < 16; i++) {
        rx |= ((d >> (2 * i)) & 1u) << i;
        ry |= ((d >> (2 * i + 1)) & 1u) << i;
    }
    *x = rx;
    *y = ry;
}

static REF_UNUSED uint64_t ref_morton2_64(uint32_t x, uint32_t y)
{
    uint64_t d = 0;
    int i;
    for (i = 0; i < 32; i++) {
        d |= ((uint64_t)(x >> i) & 1ull) << (2 * i);
        d |= ((uint64_t)(y >> i) & 1ull) << (2 * i + 1);
    }
    return d;
}

static REF_UNUSED uint64_t ref_morton3(uint32_t x, uint32_t y, uint32_t z)
{
    uint64_t d = 0;
    int i;
    x &= 0x1FFFFFu;
    y &= 0x1FFFFFu;
    z &= 0x1FFFFFu;
    for (i = 0; i < 21; i++) {
        d |= ((uint64_t)(x >> i) & 1ull) << (3 * i);
        d |= ((uint64_t)(y >> i) & 1ull) << (3 * i + 1);
        d |= ((uint64_t)(z >> i) & 1ull) << (3 * i + 2);
    }
    return d;
}

static REF_UNUSED void ref_morton3_decode(uint64_t d, uint32_t *x, uint32_t *y, uint32_t *z)
{
    uint32_t rx = 0, ry = 0, rz = 0;
    int i;
    for (i = 0; i < 21; i++) {
        rx |= (uint32_t)((d >> (3 * i)) & 1ull) << i;
        ry |= (uint32_t)((d >> (3 * i + 1)) & 1ull) << i;
        rz |= (uint32_t)((d >> (3 * i + 2)) & 1ull) << i;
    }
    *x = rx;
    *y = ry;
    *z = rz;
}

/* Wikipedia/Burkardt Hilbert: rot/d2xy/xy2d, widened to 64-bit. */
static REF_UNUSED void ref_rot(uint64_t n, uint64_t *x, uint64_t *y, uint64_t rx, uint64_t ry)
{
    if (ry == 0) {
        uint64_t t;
        if (rx == 1) { *x = n - 1 - *x; *y = n - 1 - *y; }
        t = *x; *x = *y; *y = t;
    }
}

static REF_UNUSED uint64_t ref_hilbert_xy2d(uint64_t x, uint64_t y, int bits)
{
    uint64_t d = 0, s;
    for (s = 1ull << (bits - 1); s > 0; s >>= 1) {
        uint64_t rx = (x & s) > 0;
        uint64_t ry = (y & s) > 0;
        d += s * s * ((3 * rx) ^ ry);
        ref_rot(s, &x, &y, rx, ry);
    }
    return d;
}

static REF_UNUSED void ref_hilbert_d2xy(uint64_t d, int bits, uint64_t *x, uint64_t *y)
{
    uint64_t t = d, s;
    *x = *y = 0;
    for (s = 1; s < (1ull << bits); s <<= 1) {
        uint64_t rx = 1 & (t >> 1);
        uint64_t ry = 1 & (t ^ rx);
        ref_rot(s, x, y, rx, ry);
        *x += s * rx;
        *y += s * ry;
        t >>= 2;
    }
}

/* n-D axis packing by bit-index (MSB of key first, axis 0 first in a level).
 * Independently structured from the library's nested p/i shift loop. */
static REF_UNUSED uint64_t ref_pack_axes_u64(const uint32_t *X, int n, int bits)
{
    uint64_t h = 0;
    int K = n * bits;
    int k;
    for (k = 0; k < K; k++) {
        int p = bits - 1 - (k / n);
        int i = k % n;
        uint64_t bit = ((uint64_t)X[i] >> p) & 1ull;
        h |= bit << (K - 1 - k);
    }
    return h;
}

static REF_UNUSED void ref_unpack_axes_u64(uint64_t h, int n, int bits, uint32_t *X)
{
    int K = n * bits;
    int k, i;
    for (i = 0; i < n; i++) X[i] = 0;
    for (k = 0; k < K; k++) {
        uint64_t bit = (h >> (K - 1 - k)) & 1ull;
        int p = bits - 1 - (k / n);
        int ax = k % n;
        X[ax] |= (uint32_t)bit << p;
    }
}

/* Left-aligned wide words: words[0] holds the most significant 64 key bits;
 * unused low bits of the last word are zero. */
static REF_UNUSED void ref_pack_axes_words(const uint32_t *X, int n, int bits,
                                uint64_t *words, int n_words)
{
    int K = n * bits;
    int k, w;
    for (w = 0; w < n_words; w++) words[w] = 0;
    for (k = 0; k < K; k++) {
        int p = bits - 1 - (k / n);
        int i = k % n;
        uint64_t bit = ((uint64_t)X[i] >> p) & 1ull;
        words[k / 64] |= bit << (63 - (k % 64));
    }
}

static REF_UNUSED void ref_unpack_axes_words(const uint64_t *words, int n, int bits,
                                  uint32_t *X)
{
    int K = n * bits;
    int k, i;
    for (i = 0; i < n; i++) X[i] = 0;
    for (k = 0; k < K; k++) {
        uint64_t bit = (words[k / 64] >> (63 - (k % 64))) & 1ull;
        int p = bits - 1 - (k / n);
        int ax = k % n;
        X[ax] |= (uint32_t)bit << p;
    }
}

#endif /* REF_CURVES_H */
