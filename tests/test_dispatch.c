/* Caller-owned ISA dispatch: scalar vs compiled path, allow mask, n=0/67. */
#include "curves.h"
#include "test_common.h"

#include <stdlib.h>
#include <string.h>

static void fill_xy(uint32_t *x, uint32_t *y, size_t n, uint32_t seed)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        x[i] = seed * 2654435761u + (uint32_t)i * 2246822519u;
        y[i] = seed * 3266489917u + (uint32_t)i * 40503u + 13u;
    }
}

static void test_feature_masks(void)
{
    unsigned compiled = curve_isa_compiled_features();
    unsigned runtime = curve_isa_runtime_features();
    unsigned any = curve_isa_detect(0);
    unsigned scalar = curve_isa_detect(CURVE_ISA_SCALAR);
    unsigned avx512_only = curve_isa_detect(CURVE_ISA_AVX512);
    const char *cname;

    ASSERT_TRUE(compiled & CURVE_ISA_SCALAR, "compiled includes scalar");
    ASSERT_TRUE(runtime & CURVE_ISA_SCALAR, "runtime includes scalar");
    ASSERT_TRUE(any & CURVE_ISA_SCALAR, "detect(0) includes scalar");
    ASSERT_INT_EQ(scalar, CURVE_ISA_SCALAR, "allow SCALAR is only SCALAR");
    ASSERT_TRUE(avx512_only & CURVE_ISA_SCALAR, "missing ISA still has SCALAR");
    if ((compiled & runtime & CURVE_ISA_AVX512) == 0)
        ASSERT_INT_EQ(avx512_only, CURVE_ISA_SCALAR,
                      "AVX512-only allow falls back to scalar");

    cname = curve_isa_name(any);
    ASSERT_TRUE(cname != NULL && cname[0] != '\0', "isa_name non-empty");
    ASSERT_TRUE(strcmp(curve_isa_name(0), "scalar") == 0, "name(0)=scalar");
    ASSERT_TRUE(strcmp(curve_isa_name(CURVE_ISA_NEON), "neon") == 0,
                "name(NEON)=neon");
    ASSERT_TRUE(strcmp(curve_selected_isa(NULL), "none") == 0, "NULL selected");

#if defined(CURVE_NO_SIMD)
    ASSERT_TRUE((compiled & (CURVE_ISA_NEON | CURVE_ISA_AVX2 | CURVE_ISA_AVX512)) == 0,
                "NO_SIMD compiles no vector kernels");
    ASSERT_INT_EQ(any, CURVE_ISA_SCALAR, "NO_SIMD detect is scalar");
#elif defined(CURVE_HAVE_NEON)
    ASSERT_TRUE(compiled & CURVE_ISA_NEON, "aarch64 compiles NEON");
    ASSERT_TRUE(runtime & CURVE_ISA_NEON, "aarch64 runtime NEON");
    ASSERT_TRUE(any & CURVE_ISA_NEON, "detect picks NEON");
    ASSERT_TRUE(strcmp(curve_isa_name(any), "neon") == 0, "name is neon");
#endif
}

static void test_dispatch_init_null(void)
{
    curve_dispatch_init(NULL, 0);
    ASSERT_TRUE(1, "NULL init is a no-op");
}

static int u32_eq(const uint32_t *a, const uint32_t *b, size_t n)
{
    size_t i;
    for (i = 0; i < n; i++)
        if (a[i] != b[i]) return 0;
    return 1;
}

static int u64_eq(const uint64_t *a, const uint64_t *b, size_t n)
{
    size_t i;
    for (i = 0; i < n; i++)
        if (a[i] != b[i]) return 0;
    return 1;
}

static void test_dispatch_matches_scalar(void)
{
    enum { N = 67 };
    uint32_t x[N], y[N], z[N];
    uint32_t out_s[N], out_d[N], out_p[N], rx[N], ry[N];
    uint64_t d64_s[N], d64_d[N], g[N], gout_s[N], gout_d[N];
    uint64_t m3_s[N], m3_d[N];
    curve_dispatch ds, dauto;
    size_t i;

    fill_xy(x, y, N, 7u);
    for (i = 0; i < N; i++) {
        z[i] = x[i] * 3u + 11u;
        g[i] = ((uint64_t)x[i] << 32) | y[i];
    }

    curve_dispatch_init(&ds, CURVE_ISA_SCALAR);
    curve_dispatch_init(&dauto, 0);
    ASSERT_TRUE(ds.features & CURVE_ISA_SCALAR, "scalar vtable");
    ASSERT_TRUE(strcmp(curve_selected_isa(&ds), "scalar") == 0,
                "selected scalar");
    ASSERT_TRUE(dauto.features & CURVE_ISA_SCALAR, "auto has scalar bit");

    ds.gray_encode_batch(g, gout_s, N);
    dauto.gray_encode_batch(g, gout_d, N);
    {
        uint64_t gp[N];
        curve_gray_encode_batch(g, gp, N);
        ASSERT_TRUE(u64_eq(gout_s, gp, N), "scalar gray == public");
        ASSERT_TRUE(u64_eq(gout_d, gp, N), "auto gray == public");
    }

    ds.morton2_batch(x, y, out_s, N);
    dauto.morton2_batch(x, y, out_d, N);
    curve_morton2_batch(x, y, out_p, N);
    ASSERT_TRUE(u32_eq(out_s, out_p, N), "scalar morton2 == public");
    ASSERT_TRUE(u32_eq(out_d, out_p, N), "auto morton2 == public");

    {
        uint32_t mx[N], my[N];
        size_t k;
        int ok = 1;
        ds.morton2_decode_batch(out_p, mx, my, N);
        for (k = 0; k < N; k++)
            if (mx[k] != (x[k] & 0xFFFFu) || my[k] != (y[k] & 0xFFFFu))
                ok = 0;
        ASSERT_TRUE(ok, "scalar morton2 decode");
        dauto.morton2_decode_batch(out_p, mx, my, N);
        ok = 1;
        for (k = 0; k < N; k++)
            if (mx[k] != (x[k] & 0xFFFFu) || my[k] != (y[k] & 0xFFFFu))
                ok = 0;
        ASSERT_TRUE(ok, "auto morton2 decode");
    }

    ds.morton2_64_batch(x, y, d64_s, N);
    dauto.morton2_64_batch(x, y, d64_d, N);
    {
        uint64_t p[N];
        curve_morton2_64_batch(x, y, p, N);
        ASSERT_TRUE(u64_eq(d64_s, p, N), "scalar morton2_64");
        ASSERT_TRUE(u64_eq(d64_d, p, N), "auto morton2_64");
    }

    ds.morton3_batch(x, y, z, m3_s, N);
    dauto.morton3_batch(x, y, z, m3_d, N);
    {
        uint64_t p[N];
        curve_morton3_batch(x, y, z, p, N);
        ASSERT_TRUE(u64_eq(m3_s, p, N), "scalar morton3");
        ASSERT_TRUE(u64_eq(m3_d, p, N), "auto morton3");
    }

    ds.hilbert_xy2d_batch(x, y, 8, out_s, N);
    dauto.hilbert_xy2d_batch(x, y, 8, out_d, N);
    curve_hilbert_xy2d_batch(x, y, 8, out_p, N);
    ASSERT_TRUE(u32_eq(out_s, out_p, N), "scalar hilbert xy2d");
    ASSERT_TRUE(u32_eq(out_d, out_p, N), "auto hilbert xy2d");

    ds.hilbert_d2xy_batch(out_p, 8, rx, ry, N);
    {
        uint32_t sx[N], sy[N];
        size_t k;
        int ok = 1;
        curve_hilbert_d2xy_batch(out_p, 8, sx, sy, N);
        for (k = 0; k < N; k++)
            if (rx[k] != sx[k] || ry[k] != sy[k]) ok = 0;
        ASSERT_TRUE(ok, "scalar hilbert d2xy");
        dauto.hilbert_d2xy_batch(out_p, 8, rx, ry, N);
        ok = 1;
        for (k = 0; k < N; k++)
            if (rx[k] != sx[k] || ry[k] != sy[k]) ok = 0;
        ASSERT_TRUE(ok, "auto hilbert d2xy");
    }

    ds.hilbert2_64_xy2d_batch(x, y, 20, d64_s, N);
    dauto.hilbert2_64_xy2d_batch(x, y, 20, d64_d, N);
    {
        uint64_t p[N];
        curve_hilbert2_64_xy2d_batch(x, y, 20, p, N);
        ASSERT_TRUE(u64_eq(d64_s, p, N), "scalar hilbert64");
        ASSERT_TRUE(u64_eq(d64_d, p, N), "auto hilbert64");
    }
}

static void test_dispatch_empty_and_clamp(void)
{
    curve_dispatch d;
    uint32_t x = 1, y = 2, out = 99, rx, ry;
    curve_dispatch_init(&d, CURVE_ISA_SCALAR);
    d.morton2_batch(&x, &y, &out, 0);
    ASSERT_INT_EQ(out, 99, "n=0 does not write");
    d.hilbert_xy2d_batch(&x, &y, 0, &out, 1);
    ASSERT_INT_EQ(out, curve_hilbert_xy2d(x, y, 1), "dispatch clamps bits=0");
    d.hilbert_xy2d_batch(&x, &y, 99, &out, 1);
    ASSERT_INT_EQ(out, curve_hilbert_xy2d(x, y, 16), "dispatch clamps bits=99");
    d.hilbert_d2xy_batch(&out, 8, &rx, &ry, 1);
    {
        uint32_t sx, sy;
        curve_hilbert_d2xy(out, 8, &sx, &sy);
        ASSERT_TRUE(rx == sx && ry == sy, "dispatch d2xy");
    }
}

int main(void)
{
    RUN_TEST(test_feature_masks);
    RUN_TEST(test_dispatch_init_null);
    RUN_TEST(test_dispatch_matches_scalar);
    RUN_TEST(test_dispatch_empty_and_clamp);
    TEST_SUMMARY();
}
