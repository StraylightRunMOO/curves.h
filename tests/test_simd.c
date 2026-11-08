/* SIMD Hilbert 2D batch: match scalar, tails, empty, unaligned, clamping. */
#include "curves.h"
#include "test_common.h"

#include <stdint.h>
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

static int xy2d_matches(const uint32_t *x, const uint32_t *y, int bits, size_t n)
{
    uint32_t *out;
    size_t i;
    int ok = 1;
    if (n == 0) {
        curve_hilbert_xy2d_batch(x, y, bits, NULL, 0);
        return 1;
    }
    out = (uint32_t *)malloc(n * sizeof(uint32_t));
    if (!out) return 0;
    curve_hilbert_xy2d_batch(x, y, bits, out, n);
    for (i = 0; i < n; ++i) {
        if (out[i] != curve_hilbert_xy2d(x[i], y[i], bits)) {
            ok = 0;
            break;
        }
    }
    free(out);
    return ok;
}

static int d2xy_matches(const uint32_t *d, int bits, size_t n)
{
    uint32_t *x, *y;
    size_t i;
    int ok = 1;
    if (n == 0) {
        curve_hilbert_d2xy_batch(d, bits, NULL, NULL, 0);
        return 1;
    }
    x = (uint32_t *)malloc(n * sizeof(uint32_t));
    y = (uint32_t *)malloc(n * sizeof(uint32_t));
    if (!x || !y) {
        free(x);
        free(y);
        return 0;
    }
    curve_hilbert_d2xy_batch(d, bits, x, y, n);
    for (i = 0; i < n; ++i) {
        uint32_t sx, sy;
        curve_hilbert_d2xy(d[i], bits, &sx, &sy);
        if (x[i] != sx || y[i] != sy) {
            ok = 0;
            break;
        }
    }
    free(x);
    free(y);
    return ok;
}

static int xy2d64_matches(const uint32_t *x, const uint32_t *y, int bits, size_t n)
{
    uint64_t *out;
    size_t i;
    int ok = 1;
    if (n == 0) {
        curve_hilbert2_64_xy2d_batch(x, y, bits, NULL, 0);
        return 1;
    }
    out = (uint64_t *)malloc(n * sizeof(uint64_t));
    if (!out) return 0;
    curve_hilbert2_64_xy2d_batch(x, y, bits, out, n);
    for (i = 0; i < n; ++i) {
        if (out[i] != curve_hilbert2_64_xy2d(x[i], y[i], bits)) {
            ok = 0;
            break;
        }
    }
    free(out);
    return ok;
}

static int d2xy64_matches(const uint64_t *d, int bits, size_t n)
{
    uint32_t *x, *y;
    size_t i;
    int ok = 1;
    if (n == 0) {
        curve_hilbert2_64_d2xy_batch(d, bits, NULL, NULL, 0);
        return 1;
    }
    x = (uint32_t *)malloc(n * sizeof(uint32_t));
    y = (uint32_t *)malloc(n * sizeof(uint32_t));
    if (!x || !y) {
        free(x);
        free(y);
        return 0;
    }
    curve_hilbert2_64_d2xy_batch(d, bits, x, y, n);
    for (i = 0; i < n; ++i) {
        uint32_t sx, sy;
        curve_hilbert2_64_d2xy(d[i], bits, &sx, &sy);
        if (x[i] != sx || y[i] != sy) {
            ok = 0;
            break;
        }
    }
    free(x);
    free(y);
    return ok;
}

static void test_hilbert2_batch_lengths(void)
{
    enum { MAXN = 67 };
    static const size_t ns[] = { 0, 1, 2, 3, 4, 5, 7, 8, 15, 16, 17, 31, 32, 33, 67 };
    static const int bitv[] = { 0, 1, 2, 8, 15, 16, 17, 32 };
    uint32_t x[MAXN], y[MAXN], d[MAXN];
    size_t ni, bi, i;
    fill_xy(x, y, MAXN, 1u);
    for (i = 0; i < MAXN; ++i)
        d[i] = x[i] ^ (y[i] << 1);

    ASSERT_TRUE(xy2d_matches(x, y, 16, 0), "xy2d n=0");
    ASSERT_TRUE(d2xy_matches(d, 16, 0), "d2xy n=0");

    for (ni = 0; ni < sizeof(ns) / sizeof(ns[0]); ++ni) {
        for (bi = 0; bi < sizeof(bitv) / sizeof(bitv[0]); ++bi) {
            char msg[64];
            int bits = bitv[bi];
            size_t n = ns[ni];
            snprintf(msg, sizeof(msg), "xy2d n=%zu bits=%d", n, bits);
            ASSERT_TRUE(xy2d_matches(x, y, bits, n), msg);
            snprintf(msg, sizeof(msg), "d2xy n=%zu bits=%d", n, bits);
            ASSERT_TRUE(d2xy_matches(d, bits, n), msg);
        }
    }
}

static void test_hilbert2_64_batch_lengths(void)
{
    enum { MAXN = 67 };
    static const size_t ns[] = { 0, 1, 2, 3, 4, 5, 7, 8, 15, 16, 17, 67 };
    static const int bitv[] = { 0, 1, 16, 17, 31, 32, 33 };
    uint32_t x[MAXN], y[MAXN];
    uint64_t d[MAXN];
    size_t ni, bi, i;
    fill_xy(x, y, MAXN, 9u);
    for (i = 0; i < MAXN; ++i)
        d[i] = ((uint64_t)x[i] << 32) | y[i];

    ASSERT_TRUE(xy2d64_matches(x, y, 32, 0), "xy2d64 n=0");
    ASSERT_TRUE(d2xy64_matches(d, 32, 0), "d2xy64 n=0");

    for (ni = 0; ni < sizeof(ns) / sizeof(ns[0]); ++ni) {
        for (bi = 0; bi < sizeof(bitv) / sizeof(bitv[0]); ++bi) {
            char msg[64];
            int bits = bitv[bi];
            size_t n = ns[ni];
            snprintf(msg, sizeof(msg), "xy2d64 n=%zu bits=%d", n, bits);
            ASSERT_TRUE(xy2d64_matches(x, y, bits, n), msg);
            snprintf(msg, sizeof(msg), "d2xy64 n=%zu bits=%d", n, bits);
            ASSERT_TRUE(d2xy64_matches(d, bits, n), msg);
        }
    }
}

static void test_hilbert2_batch_unaligned(void)
{
    enum { N = 19 };
    uint32_t bufx[N + 1], bufy[N + 1], bufo[N + 1], bufd[N + 1];
    uint32_t *x = bufx + 1, *y = bufy + 1, *o = bufo + 1, *d = bufd + 1;
    uint32_t rx[N], ry[N];
    size_t i;
    int ok;

    fill_xy(x, y, N, 42u);
    curve_hilbert_xy2d_batch(x, y, 12, o, N);
    ok = 1;
    for (i = 0; i < N; ++i)
        if (o[i] != curve_hilbert_xy2d(x[i], y[i], 12)) ok = 0;
    ASSERT_TRUE(ok, "unaligned xy2d +1");

    memcpy(d, o, N * sizeof(uint32_t));
    curve_hilbert_d2xy_batch(d, 12, rx, ry, N);
    ok = 1;
    for (i = 0; i < N; ++i) {
        uint32_t sx, sy;
        curve_hilbert_d2xy(d[i], 12, &sx, &sy);
        if (rx[i] != sx || ry[i] != sy) ok = 0;
    }
    ASSERT_TRUE(ok, "unaligned d2xy (out aligned, in +1)");
}

static void test_hilbert2_batch_roundtrip(void)
{
    enum { N = 64 };
    uint32_t x[N], y[N], d[N], rx[N], ry[N];
    uint32_t mask;
    size_t i;
    int bits, ok;

    for (bits = 1; bits <= 16; ++bits) {
        mask = (bits >= 32) ? 0xFFFFFFFFu : ((1u << bits) - 1u);
        fill_xy(x, y, N, (uint32_t)bits * 17u);
        for (i = 0; i < N; ++i) {
            x[i] &= mask;
            y[i] &= mask;
        }
        curve_hilbert_xy2d_batch(x, y, bits, d, N);
        curve_hilbert_d2xy_batch(d, bits, rx, ry, N);
        ok = 1;
        for (i = 0; i < N; ++i)
            if (rx[i] != x[i] || ry[i] != y[i]) ok = 0;
        ASSERT_TRUE(ok, "hilbert2 batch roundtrip");
    }
}

static void test_hilbert2_64_batch_roundtrip(void)
{
    enum { N = 32 };
    uint32_t x[N], y[N], rx[N], ry[N];
    uint64_t d[N];
    uint32_t mask;
    size_t i;
    int bits, ok;

    for (bits = 1; bits <= 32; bits += 5) {
        mask = (bits >= 32) ? 0xFFFFFFFFu : ((1u << bits) - 1u);
        fill_xy(x, y, N, (uint32_t)bits * 91u);
        for (i = 0; i < N; ++i) {
            x[i] &= mask;
            y[i] &= mask;
        }
        curve_hilbert2_64_xy2d_batch(x, y, bits, d, N);
        curve_hilbert2_64_d2xy_batch(d, bits, rx, ry, N);
        ok = 1;
        for (i = 0; i < N; ++i)
            if (rx[i] != x[i] || ry[i] != y[i]) ok = 0;
        ASSERT_TRUE(ok, "hilbert2_64 batch roundtrip");
    }
}

static void test_compiled_isa(void)
{
    const char *a = curve_compiled_isa();
    const char *b = curve_simd_isa();
    ASSERT_TRUE(a != NULL && a[0] != '\0', "compiled_isa non-empty");
    ASSERT_TRUE(b != NULL && strcmp(a, b) == 0, "simd_isa aliases compiled_isa");
#if defined(CURVE_NO_SIMD)
    ASSERT_TRUE(strcmp(a, "scalar") == 0 || strcmp(a, "bmi2") == 0,
                "CURVE_NO_SIMD is scalar or bmi2");
#endif
}

int main(void)
{
    RUN_TEST(test_hilbert2_batch_lengths);
    RUN_TEST(test_hilbert2_64_batch_lengths);
    RUN_TEST(test_hilbert2_batch_unaligned);
    RUN_TEST(test_hilbert2_batch_roundtrip);
    RUN_TEST(test_hilbert2_64_batch_roundtrip);
    RUN_TEST(test_compiled_isa);
    TEST_SUMMARY();
}
