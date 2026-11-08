#include "curves.h"
#include "test_common.h"

#include <string.h>

static double dnan(void)
{
    volatile double z = 0.0;
    return z / z;
}

static double dpinf(void)
{
    volatile double z = 0.0;
    return 1.0 / z;
}

static double dninf(void)
{
    volatile double z = 0.0;
    return -1.0 / z;
}

static float fnan(void)
{
    volatile float z = 0.0f;
    return z / z;
}

static float fpinf(void)
{
    volatile float z = 0.0f;
    return 1.0f / z;
}

static int keys_equal(const uint64_t *a, const uint64_t *b, size_t n)
{
    size_t i;
    for (i = 0; i < n; i++)
        if (a[i] != b[i]) return 0;
    return 1;
}

static void test_known_floor(void)
{
    uint32_t q = 99u;
    /* [0,4], bits=2, n_bins=4, width=1: [0,1) [1,2) [2,3) [3,4]. */
    ASSERT_INT_EQ(curve_quantize1_f64(0.0, 0.0, 4.0, 2, CURVE_QUANT_FLOOR, &q),
                  CURVE_OK, "floor 0");
    ASSERT_INT_EQ(q, 0u, "0 -> 0");
    ASSERT_INT_EQ(curve_quantize1_f64(0.999, 0.0, 4.0, 2, CURVE_QUANT_FLOOR, &q),
                  CURVE_OK, "floor 0.999");
    ASSERT_INT_EQ(q, 0u, "0.999 -> 0");
    ASSERT_INT_EQ(curve_quantize1_f64(1.0, 0.0, 4.0, 2, CURVE_QUANT_FLOOR, &q),
                  CURVE_OK, "floor 1");
    ASSERT_INT_EQ(q, 1u, "1 -> 1");
    ASSERT_INT_EQ(curve_quantize1_f64(3.999, 0.0, 4.0, 2, CURVE_QUANT_FLOOR, &q),
                  CURVE_OK, "floor 3.999");
    ASSERT_INT_EQ(q, 3u, "3.999 -> 3");
    ASSERT_INT_EQ(curve_quantize1_f64(4.0, 0.0, 4.0, 2, CURVE_QUANT_FLOOR, &q),
                  CURVE_OK, "floor 4");
    ASSERT_INT_EQ(q, 3u, "max -> last");
}

static void test_known_nearest(void)
{
    uint32_t q = 99u;
    /* [0,1], bits=2, last=3, grid 0, 1/3, 2/3, 1. */
    ASSERT_INT_EQ(curve_quantize1_f64(0.0, 0.0, 1.0, 2, CURVE_QUANT_NEAREST, &q),
                  CURVE_OK, "near 0");
    ASSERT_INT_EQ(q, 0u, "0 -> 0");
    ASSERT_INT_EQ(curve_quantize1_f64(1.0, 0.0, 1.0, 2, CURVE_QUANT_NEAREST, &q),
                  CURVE_OK, "near 1");
    ASSERT_INT_EQ(q, 3u, "1 -> 3");
    /* 0.5 is equidistant between 1/3 and 2/3; tie toward +inf -> 2. */
    ASSERT_INT_EQ(curve_quantize1_f64(0.5, 0.0, 1.0, 2, CURVE_QUANT_NEAREST, &q),
                  CURVE_OK, "near 0.5 tie");
    ASSERT_INT_EQ(q, 2u, "0.5 tie -> 2");
    /* [0,7] bits=3 is the integer grid; 3.5 ties 3 vs 4 -> 4. */
    ASSERT_INT_EQ(curve_quantize1_f64(3.5, 0.0, 7.0, 3, CURVE_QUANT_NEAREST, &q),
                  CURVE_OK, "near 3.5");
    ASSERT_INT_EQ(q, 4u, "3.5 tie -> 4");
    ASSERT_INT_EQ(curve_quantize1_f64(3.0, 0.0, 7.0, 3, CURVE_QUANT_NEAREST, &q),
                  CURVE_OK, "near 3");
    ASSERT_INT_EQ(q, 3u, "3 -> 3");
}

static void test_nan_inf_zero(void)
{
    uint32_t q = 99u;
    size_t n_lo = 99, n_hi = 99, n_nan = 99;
    double v[6], mn[6], mx[6];
    uint32_t out[6];
    size_t i;
    double mz = 0.0, nz;

    nz = -mz; /* -0.0 */
    ASSERT_INT_EQ(curve_quantize1_f64(dnan(), 0.0, 1.0, 4, CURVE_QUANT_FLOOR, &q),
                  CURVE_OK, "nan ok");
    ASSERT_INT_EQ(q, 0u, "nan -> 0");
    ASSERT_INT_EQ(curve_quantize1_f64(dninf(), 0.0, 1.0, 4, CURVE_QUANT_FLOOR, &q),
                  CURVE_OK, "-inf ok");
    ASSERT_INT_EQ(q, 0u, "-inf -> 0");
    ASSERT_INT_EQ(curve_quantize1_f64(dpinf(), 0.0, 1.0, 4, CURVE_QUANT_FLOOR, &q),
                  CURVE_OK, "+inf ok");
    ASSERT_INT_EQ(q, 15u, "+inf -> last");
    ASSERT_INT_EQ(curve_quantize1_f64(nz, 0.0, 1.0, 4, CURVE_QUANT_NEAREST, &q),
                  CURVE_OK, "-0 ok");
    ASSERT_INT_EQ(q, 0u, "-0 -> 0");

    for (i = 0; i < 6; i++) {
        mn[i] = 0.0;
        mx[i] = 1.0;
    }
    v[0] = dnan();
    v[1] = dninf();
    v[2] = -0.5;
    v[3] = 0.5;
    v[4] = 2.0;
    v[5] = dpinf();
    ASSERT_INT_EQ(curve_quantize_f64(v, 6, mn, mx, 4, CURVE_QUANT_FLOOR,
                                     out, 6, &n_lo, &n_hi, &n_nan),
                  CURVE_OK, "batch specials");
    ASSERT_INT_EQ(n_nan, 1u, "one nan");
    ASSERT_INT_EQ(n_lo, 2u, "-inf and -0.5");
    ASSERT_INT_EQ(n_hi, 2u, "2.0 and +inf");
    ASSERT_INT_EQ(out[0], 0u, "nan bin");
    ASSERT_INT_EQ(out[1], 0u, "-inf bin");
    ASSERT_INT_EQ(out[2], 0u, "below bin");
    ASSERT_INT_EQ(out[4], 15u, "above bin");
    ASSERT_INT_EQ(out[5], 15u, "+inf bin");

    n_lo = 77;
    n_hi = 77;
    n_nan = 77;
    ASSERT_INT_EQ(curve_quantize_f64(v, 6, mn, mx, 4, CURVE_QUANT_FLOOR,
                                     out, 6, &n_lo, &n_hi, &n_nan),
                  CURVE_OK, "counters set not accumulated");
    ASSERT_INT_EQ(n_lo, 2u, "lo reset");
    ASSERT_INT_EQ(n_hi, 2u, "hi reset");
    ASSERT_INT_EQ(n_nan, 1u, "nan reset");

    ASSERT_INT_EQ(curve_quantize_f64(v, 6, mn, mx, 4, CURVE_QUANT_FLOOR,
                                     out, 6, NULL, NULL, NULL),
                  CURVE_OK, "null counters");
}

static void test_zero_width_and_errors(void)
{
    uint32_t q = 99u;
    size_t n_lo = 0, n_hi = 0, n_nan = 0;
    double v[3] = { 1.0, 2.0, 3.0 };
    double mn[3] = { 5.0, 5.0, 5.0 };
    double mx[3] = { 5.0, 5.0, 5.0 };
    uint32_t out[3];
    double bad_mn[1] = { 1.0 }, bad_mx[1] = { 0.0 };
    double nan_mn[1], inf_mx[1];
    uint32_t sentinel = 0xABCDu;

    ASSERT_INT_EQ(curve_quantize1_f64(5.0, 5.0, 5.0, 8, CURVE_QUANT_FLOOR, &q),
                  CURVE_OK, "zero width");
    ASSERT_INT_EQ(q, 0u, "zero width -> 0");
    ASSERT_INT_EQ(curve_quantize1_f64(dpinf(), 5.0, 5.0, 8, CURVE_QUANT_NEAREST,
                                      &q),
                  CURVE_OK, "zero width +inf");
    ASSERT_INT_EQ(q, 0u, "zero width +inf -> 0");

    ASSERT_INT_EQ(curve_quantize_f64(v, 3, mn, mx, 4, CURVE_QUANT_FLOOR,
                                     out, 3, &n_lo, &n_hi, &n_nan),
                  CURVE_OK, "zero width batch");
    ASSERT_INT_EQ(out[0], 0u, "zw 0");
    ASSERT_INT_EQ(out[1], 0u, "zw 1");
    ASSERT_INT_EQ(n_lo, 3u, "all below the point");
    ASSERT_INT_EQ(n_hi, 0u, "none above");

    q = sentinel;
    ASSERT_INT_EQ(curve_quantize1_f64(0.5, 1.0, 0.0, 4, CURVE_QUANT_FLOOR, &q),
                  CURVE_EINVAL, "min>max");
    ASSERT_INT_EQ(q, sentinel, "no write on EINVAL");
    nan_mn[0] = dnan();
    inf_mx[0] = 1.0;
    ASSERT_INT_EQ(curve_quantize_f64(v, 1, nan_mn, inf_mx, 4,
                                     CURVE_QUANT_FLOOR, out, 1,
                                     NULL, NULL, NULL),
                  CURVE_EINVAL, "nan bound");
    inf_mx[0] = dpinf();
    nan_mn[0] = 0.0;
    ASSERT_INT_EQ(curve_quantize_f64(v, 1, nan_mn, inf_mx, 4,
                                     CURVE_QUANT_FLOOR, out, 1,
                                     NULL, NULL, NULL),
                  CURVE_EINVAL, "inf bound");
    ASSERT_INT_EQ(curve_quantize_f64(v, 1, bad_mn, bad_mx, 4,
                                     CURVE_QUANT_FLOOR, out, 1,
                                     NULL, NULL, NULL),
                  CURVE_EINVAL, "hi<lo array");
    ASSERT_INT_EQ(curve_quantize1_f64(0.5, 0.0, 1.0, 0, CURVE_QUANT_FLOOR, &q),
                  CURVE_EINVAL, "bits 0");
    ASSERT_INT_EQ(curve_quantize1_f64(0.5, 0.0, 1.0, 33, CURVE_QUANT_FLOOR, &q),
                  CURVE_EINVAL, "bits 33");
    ASSERT_INT_EQ(curve_quantize1_f64(0.5, 0.0, 1.0, 4,
                                      (curve_quant_round)2, &q),
                  CURVE_EINVAL, "bad rnd");
    ASSERT_INT_EQ(curve_quantize1_f64(0.5, 0.0, 1.0, 4, CURVE_QUANT_FLOOR, NULL),
                  CURVE_EINVAL, "null out");
    ASSERT_INT_EQ(curve_quantize_f64(v, 3, mn, mx, 4, CURVE_QUANT_FLOOR,
                                     out, 2, NULL, NULL, NULL),
                  CURVE_ENOSPACE, "short cap");
    ASSERT_INT_EQ(curve_quantize_f64(NULL, 0, NULL, NULL, 4, CURVE_QUANT_FLOOR,
                                     NULL, 0, &n_lo, &n_hi, &n_nan),
                  CURVE_OK, "n=0");
    ASSERT_INT_EQ(n_lo, 0u, "n=0 lo");
    ASSERT_INT_EQ(n_hi, 0u, "n=0 hi");
    ASSERT_INT_EQ(n_nan, 0u, "n=0 nan");
}

static void test_bits32_and_monotonic(void)
{
    uint32_t q = 0, prev = 0;
    int i, ok = 1;
    double x;

    ASSERT_INT_EQ(curve_quantize1_f64(0.0, 0.0, 1.0, 32, CURVE_QUANT_FLOOR, &q),
                  CURVE_OK, "b32 floor 0");
    ASSERT_INT_EQ(q, 0u, "b32 0");
    ASSERT_INT_EQ(curve_quantize1_f64(1.0, 0.0, 1.0, 32, CURVE_QUANT_FLOOR, &q),
                  CURVE_OK, "b32 floor 1");
    ASSERT_INT_EQ(q, 0xFFFFFFFFu, "b32 floor max");
    ASSERT_INT_EQ(curve_quantize1_f64(1.0, 0.0, 1.0, 32, CURVE_QUANT_NEAREST,
                                      &q),
                  CURVE_OK, "b32 near 1");
    ASSERT_INT_EQ(q, 0xFFFFFFFFu, "b32 near max");
    ASSERT_INT_EQ(curve_quantize1_f64(0.5, 0.0, 1.0, 32, CURVE_QUANT_NEAREST,
                                      &q),
                  CURVE_OK, "b32 near mid");
    ASSERT_INT_EQ(q, 0x80000000u, "b32 mid -> 2^31");

    prev = 0;
    for (i = 0; i <= 200; i++) {
        x = (double)i / 200.0 * 3.0 - 0.5; /* [-0.5, 2.5] over [0,2] */
        ASSERT_INT_EQ(curve_quantize1_f64(x, 0.0, 2.0, 8, CURVE_QUANT_FLOOR, &q),
                      CURVE_OK, "mono floor");
        if (i > 0 && q < prev) ok = 0;
        prev = q;
    }
    ASSERT_TRUE(ok, "floor monotonic");

    ok = 1;
    prev = 0;
    for (i = 0; i <= 200; i++) {
        x = (double)i / 200.0;
        ASSERT_INT_EQ(curve_quantize1_f64(x, 0.0, 1.0, 12, CURVE_QUANT_NEAREST,
                                          &q),
                      CURVE_OK, "mono near");
        if (i > 0 && q < prev) ok = 0;
        prev = q;
    }
    ASSERT_TRUE(ok, "nearest monotonic");
}

static void test_negative_range(void)
{
    uint32_t q = 0;
    /* [-2,2] bits=2 floor: width=1, bins [-2,-1) [-1,0) [0,1) [1,2]. */
    ASSERT_INT_EQ(curve_quantize1_f64(-2.0, -2.0, 2.0, 2, CURVE_QUANT_FLOOR, &q),
                  CURVE_OK, "neg lo");
    ASSERT_INT_EQ(q, 0u, "-2 -> 0");
    ASSERT_INT_EQ(curve_quantize1_f64(-1.0, -2.0, 2.0, 2, CURVE_QUANT_FLOOR, &q),
                  CURVE_OK, "neg -1");
    ASSERT_INT_EQ(q, 1u, "-1 -> 1");
    ASSERT_INT_EQ(curve_quantize1_f64(0.0, -2.0, 2.0, 2, CURVE_QUANT_FLOOR, &q),
                  CURVE_OK, "neg 0");
    ASSERT_INT_EQ(q, 2u, "0 -> 2");
    ASSERT_INT_EQ(curve_quantize1_f64(2.0, -2.0, 2.0, 2, CURVE_QUANT_FLOOR, &q),
                  CURVE_OK, "neg hi");
    ASSERT_INT_EQ(q, 3u, "2 -> 3");
}

static void test_dequantize(void)
{
    uint32_t q;
    double y = 99.0;
    float yf = 99.0f;
    uint32_t qs[4] = { 0, 1, 2, 3 };
    double mn[4] = { 0, 0, 0, 0 }, mx[4] = { 4, 4, 4, 4 }, out[4];

    ASSERT_INT_EQ(curve_quantize1_f64(0.0, 0.0, 1.0, 4, CURVE_QUANT_NEAREST, &q),
                  CURVE_OK, "dq q lo");
    ASSERT_INT_EQ(curve_dequantize1_f64(q, 0.0, 1.0, 4, CURVE_QUANT_NEAREST, &y),
                  CURVE_OK, "dq lo");
    ASSERT_TRUE(y == 0.0, "nearest reconstructs min");
    ASSERT_INT_EQ(curve_quantize1_f64(1.0, 0.0, 1.0, 4, CURVE_QUANT_NEAREST, &q),
                  CURVE_OK, "dq q hi");
    ASSERT_INT_EQ(curve_dequantize1_f64(q, 0.0, 1.0, 4, CURVE_QUANT_NEAREST, &y),
                  CURVE_OK, "dq hi");
    ASSERT_TRUE(y == 1.0, "nearest reconstructs max");

    ASSERT_INT_EQ(curve_dequantize1_f64(0, 0.0, 4.0, 2, CURVE_QUANT_FLOOR, &y),
                  CURVE_OK, "floor center 0");
    ASSERT_TRUE(y == 0.5, "bin0 center");
    ASSERT_INT_EQ(curve_dequantize1_f64(3, 0.0, 4.0, 2, CURVE_QUANT_FLOOR, &y),
                  CURVE_OK, "floor center last");
    ASSERT_TRUE(y == 3.5, "bin3 center");

    ASSERT_INT_EQ(curve_dequantize_f64(qs, 4, mn, mx, 2, CURVE_QUANT_FLOOR,
                                       out, 4),
                  CURVE_OK, "dequant array");
    ASSERT_TRUE(out[0] == 0.5 && out[3] == 3.5, "array centers");

    ASSERT_INT_EQ(curve_dequantize1_f64(4, 0.0, 4.0, 2, CURVE_QUANT_FLOOR, &y),
                  CURVE_EINVAL, "q > last");
    ASSERT_INT_EQ(curve_dequantize1_f32(0, 0.0f, 1.0f, 4, CURVE_QUANT_NEAREST,
                                        &yf),
                  CURVE_OK, "f32 dequant");
    ASSERT_TRUE(yf == 0.0f, "f32 min");
}

static void test_f32_matches_f64(void)
{
    float xf[5] = { 0.0f, 0.25f, 0.5f, 0.75f, 1.0f };
    float mn[5], mx[5];
    double xd[5], dmn[5], dmx[5];
    uint32_t qf[5], qd[5];
    size_t i;
    for (i = 0; i < 5; i++) {
        mn[i] = 0.0f;
        mx[i] = 1.0f;
        xd[i] = (double)xf[i];
        dmn[i] = 0.0;
        dmx[i] = 1.0;
    }
    ASSERT_INT_EQ(curve_quantize_f32(xf, 5, mn, mx, 8, CURVE_QUANT_NEAREST,
                                     qf, 5, NULL, NULL, NULL),
                  CURVE_OK, "f32");
    ASSERT_INT_EQ(curve_quantize_f64(xd, 5, dmn, dmx, 8, CURVE_QUANT_NEAREST,
                                     qd, 5, NULL, NULL, NULL),
                  CURVE_OK, "f64");
    for (i = 0; i < 5; i++)
        ASSERT_INT_EQ(qf[i], qd[i], "f32 matches promoted f64");

    {
        uint32_t q = 99u;
        ASSERT_INT_EQ(curve_quantize1_f32(fnan(), 0.0f, 1.0f, 4,
                                          CURVE_QUANT_FLOOR, &q),
                      CURVE_OK, "f32 nan");
        ASSERT_INT_EQ(q, 0u, "f32 nan -> 0");
        ASSERT_INT_EQ(curve_quantize1_f32(fpinf(), 0.0f, 1.0f, 4,
                                          CURVE_QUANT_FLOOR, &q),
                      CURVE_OK, "f32 +inf");
        ASSERT_INT_EQ(q, 15u, "f32 +inf -> last");
    }
}

static void test_fused_equals_unfused(void)
{
    enum { N = 256, B = 4 };
    double v[N], mn[N], mx[N];
    float vf[N], mnf[N], mxf[N];
    uint32_t q[N], scratch[N];
    uint64_t km[16], kh[16], kfm[16], kfh[16], k32m[16], k32h[16];
    size_t i, W;
    size_t n_lo = 0, n_hi = 0, n_nan = 0;

    W = curve_key_words_required(N, B);
    ASSERT_INT_EQ(W, 16u, "256x4 words");
    for (i = 0; i < (size_t)N; i++) {
        v[i] = (double)i * 0.013 + 0.02;
        mn[i] = 0.0;
        mx[i] = 4.0;
        vf[i] = (float)v[i];
        mnf[i] = 0.0f;
        mxf[i] = 4.0f;
    }
    v[0] = -1.0;
    v[1] = 9.0;
    v[2] = dnan();
    vf[0] = -1.0f;
    vf[1] = 9.0f;
    vf[2] = fnan();

    ASSERT_INT_EQ(curve_quantize_f64(v, N, mn, mx, B, CURVE_QUANT_NEAREST,
                                     q, N, &n_lo, &n_hi, &n_nan),
                  CURVE_OK, "quantize 256");
    ASSERT_TRUE(n_lo >= 1 && n_hi >= 1 && n_nan >= 1, "saturations");
    ASSERT_INT_EQ(curve_mortonnd_encode_words(q, N, B, km, 16),
                  CURVE_OK, "morton unfused");
    ASSERT_INT_EQ(curve_hilbertnd_encode_words(q, N, B, kh, 16, scratch, N),
                  CURVE_OK, "hilbert unfused");

    ASSERT_INT_EQ(curve_mortonnd_encode_f64(v, N, mn, mx, B, CURVE_QUANT_NEAREST,
                                            kfm, 16, scratch, N,
                                            NULL, NULL, NULL),
                  CURVE_OK, "morton fused");
    ASSERT_TRUE(keys_equal(km, kfm, W), "morton fused==unfused");

    ASSERT_INT_EQ(curve_hilbertnd_encode_f64(v, N, mn, mx, B,
                                             CURVE_QUANT_NEAREST,
                                             kfh, 16, scratch, N,
                                             NULL, NULL, NULL),
                  CURVE_OK, "hilbert fused");
    ASSERT_TRUE(keys_equal(kh, kfh, W), "hilbert fused==unfused");

    ASSERT_INT_EQ(curve_mortonnd_encode_f32(vf, N, mnf, mxf, B,
                                            CURVE_QUANT_NEAREST,
                                            k32m, 16, scratch, N,
                                            NULL, NULL, NULL),
                  CURVE_OK, "morton f32 fused");
    ASSERT_INT_EQ(curve_hilbertnd_encode_f32(vf, N, mnf, mxf, B,
                                             CURVE_QUANT_NEAREST,
                                             k32h, 16, scratch, N,
                                             NULL, NULL, NULL),
                  CURVE_OK, "hilbert f32 fused");
}

static void test_aos_batch(void)
{
    enum { NV = 5, D = 8, B = 4 };
    double vals[NV * D], mn[D], mx[D];
    uint64_t keys[NV], loop[NV];
    uint32_t scratch[D], q[D];
    size_t i, d, n_lo = 0, n_hi = 0, n_nan = 0;
    size_t l_lo = 0, l_hi = 0, l_nan = 0;

    for (d = 0; d < (size_t)D; d++) {
        mn[d] = 0.0;
        mx[d] = 1.0;
    }
    for (i = 0; i < (size_t)NV; i++) {
        for (d = 0; d < (size_t)D; d++)
            vals[i * D + d] = (double)((i + d) % 11) / 10.0;
    }
    vals[2] = -0.25;
    vals[D + 1] = 2.0;

    ASSERT_INT_EQ(curve_mortonnd_encode_f64_aos(vals, NV, D, D, mn, mx, B,
                                                CURVE_QUANT_FLOOR,
                                                keys, 1, scratch, D,
                                                &n_lo, &n_hi, &n_nan),
                  CURVE_OK, "aos morton");
    for (i = 0; i < (size_t)NV; i++) {
        size_t c_lo = 0, c_hi = 0, c_nan = 0;
        ASSERT_INT_EQ(curve_quantize_f64(vals + i * D, D, mn, mx, B,
                                         CURVE_QUANT_FLOOR, q, D,
                                         &c_lo, &c_hi, &c_nan),
                      CURVE_OK, "loop q");
        ASSERT_INT_EQ(curve_mortonnd_encode_words(q, D, B, &loop[i], 1),
                      CURVE_OK, "loop enc");
        l_lo += c_lo;
        l_hi += c_hi;
        l_nan += c_nan;
    }
    ASSERT_TRUE(keys_equal(keys, loop, NV), "aos==loop morton");
    ASSERT_INT_EQ(n_lo, l_lo, "aos lo");
    ASSERT_INT_EQ(n_hi, l_hi, "aos hi");
    ASSERT_INT_EQ(n_nan, l_nan, "aos nan");

    ASSERT_INT_EQ(curve_hilbertnd_encode_f64_aos(vals, NV, D, D, mn, mx, B,
                                                 CURVE_QUANT_FLOOR,
                                                 keys, 1, scratch, D,
                                                 NULL, NULL, NULL),
                  CURVE_OK, "aos hilbert");
    for (i = 0; i < (size_t)NV; i++) {
        ASSERT_INT_EQ(curve_quantize_f64(vals + i * D, D, mn, mx, B,
                                         CURVE_QUANT_FLOOR, q, D,
                                         NULL, NULL, NULL),
                      CURVE_OK, "h loop q");
        ASSERT_INT_EQ(curve_hilbertnd_encode_words(q, D, B, &loop[i], 1,
                                                   scratch, D),
                      CURVE_OK, "h loop enc");
    }
    ASSERT_TRUE(keys_equal(keys, loop, NV), "aos==loop hilbert");

    ASSERT_INT_EQ(curve_mortonnd_encode_f64_aos(NULL, 0, D, D, NULL, NULL, B,
                                                CURVE_QUANT_FLOOR,
                                                NULL, 1, NULL, 0,
                                                &n_lo, &n_hi, &n_nan),
                  CURVE_OK, "aos n_vecs=0");
    ASSERT_INT_EQ(n_lo + n_hi + n_nan, 0u, "aos empty counts");
    ASSERT_INT_EQ(curve_mortonnd_encode_f64_aos(vals, NV, D - 1, D, mn, mx, B,
                                                CURVE_QUANT_FLOOR,
                                                keys, 1, scratch, D,
                                                NULL, NULL, NULL),
                  CURVE_EINVAL, "short value stride");
}

static void test_fused_errors(void)
{
    double v[2] = { 0.5, 0.5 }, mn[2] = { 0, 0 }, mx[2] = { 1, 1 };
    uint64_t key[1];
    uint32_t scratch[2];
    uint64_t sentinel = 0xDEADBEEFCAFEBABEull;

    key[0] = sentinel;
    ASSERT_INT_EQ(curve_mortonnd_encode_f64(v, 2, mn, mx, 8, CURVE_QUANT_FLOOR,
                                            key, 0, scratch, 2,
                                            NULL, NULL, NULL),
                  CURVE_ENOSPACE, "short key");
    ASSERT_TRUE(key[0] == sentinel, "key unchanged");
    ASSERT_INT_EQ(curve_mortonnd_encode_f64(v, 2, mn, mx, 8, CURVE_QUANT_FLOOR,
                                            key, 1, scratch, 1,
                                            NULL, NULL, NULL),
                  CURVE_ENOSPACE, "short scratch");
    mx[0] = -1.0;
    ASSERT_INT_EQ(curve_hilbertnd_encode_f64(v, 2, mn, mx, 8, CURVE_QUANT_FLOOR,
                                             key, 1, scratch, 2,
                                             NULL, NULL, NULL),
                  CURVE_EINVAL, "bad bounds fused");
    ASSERT_INT_EQ(curve_mortonnd_encode_f64(NULL, 2, mn, mx, 8,
                                            CURVE_QUANT_FLOOR,
                                            key, 1, scratch, 2,
                                            NULL, NULL, NULL),
                  CURVE_EINVAL, "null values");
}

int main(void)
{
    RUN_TEST(test_known_floor);
    RUN_TEST(test_known_nearest);
    RUN_TEST(test_nan_inf_zero);
    RUN_TEST(test_zero_width_and_errors);
    RUN_TEST(test_bits32_and_monotonic);
    RUN_TEST(test_negative_range);
    RUN_TEST(test_dequantize);
    RUN_TEST(test_f32_matches_f64);
    RUN_TEST(test_fused_equals_unfused);
    RUN_TEST(test_aos_batch);
    RUN_TEST(test_fused_errors);
    TEST_SUMMARY();
}
