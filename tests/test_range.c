#include "curves.h"
#include "test_common.h"

#include <stdlib.h>
#include <string.h>

static size_t g_sort_words;

static int cmp_keys(const void *pa, const void *pb)
{
    return curve_key_compare((const uint64_t *)pa, (const uint64_t *)pb,
                             g_sort_words);
}

static uint32_t rng_u32(uint32_t *s)
{
    *s = *s * 1664525u + 1013904223u;
    return *s;
}

static int encode_pt(int hilbert, size_t n, unsigned bits,
                     const uint32_t *coords, uint64_t *key)
{
    uint32_t scratch[8];
    if (hilbert)
        return curve_hilbertnd_encode_words(coords, n, bits, key, 1,
                                            scratch, 8) == CURVE_OK;
    return curve_mortonnd_encode_words(coords, n, bits, key, 1) == CURVE_OK;
}

static int point_in_box(const uint32_t *c, const uint32_t *bmin,
                        const uint32_t *bmax, size_t n)
{
    size_t i;
    for (i = 0; i < n; i++)
        if (c[i] < bmin[i] || c[i] > bmax[i]) return 0;
    return 1;
}

static int key_in_span(const uint64_t *key, size_t W,
                       const uint64_t *lo, const uint64_t *hi, int unb)
{
    if (curve_key_compare(key, lo, W) < 0) return 0;
    if (unb) return 1;
    return curve_key_compare(key, hi, W) < 0;
}

/* Exact oracle: encode every in-box point, sort, merge successor runs. */
static int oracle_intervals(int hilbert, size_t n, unsigned bits,
                            const uint32_t *bmin, const uint32_t *bmax,
                            uint64_t *lo, uint64_t *hi, int *unb,
                            size_t cap, size_t *out_n)
{
    uint32_t coords[8], t, side, total, v;
    uint64_t keys[4096], succ[1];
    size_t i, nk = 0, nint = 0, d, W, K;
    curve_status st;

    if (n > 8u) return 0;
    side = 1u << bits;
    total = 1;
    for (d = 0; d < n; d++) total *= side;
    if (total > 4096u) return 0;
    K = n * (size_t)bits;
    W = (K + 63u) / 64u;
    for (t = 0; t < total; t++) {
        v = t;
        for (d = 0; d < n; d++) {
            coords[d] = v % side;
            v /= side;
        }
        if (!point_in_box(coords, bmin, bmax, n)) continue;
        if (!encode_pt(hilbert, n, bits, coords, &keys[nk])) return 0;
        nk++;
    }
    g_sort_words = W;
    qsort(keys, nk, sizeof(uint64_t), cmp_keys);
    i = 0;
    while (i < nk) {
        size_t j = i;
        while (j + 1u < nk) {
            st = curve_key_prefix_successor(&keys[j], K, K, succ, 1);
            if (st == CURVE_ERANGE) break;
            if (st != CURVE_OK) return 0;
            if (succ[0] != keys[j + 1u]) break;
            j++;
        }
        if (nint < cap) {
            lo[nint] = keys[i];
            st = curve_key_prefix_successor(&keys[j], K, K, succ, 1);
            if (st == CURVE_ERANGE) {
                hi[nint] = 0;
                unb[nint] = 1;
            } else if (st == CURVE_OK) {
                hi[nint] = succ[0];
                unb[nint] = 0;
            } else {
                return 0;
            }
        }
        nint++;
        i = j + 1u;
    }
    *out_n = nint;
    return 1;
}

static int check_box(int hilbert, size_t n, unsigned bits,
                     const uint32_t *bmin, const uint32_t *bmax)
{
    uint64_t scratch[128], out_lo[64], out_hi[64], olo[64], ohi[64];
    int out_unb[64], ounb[64];
    uint32_t axes[40];
    size_t count = 0, ocount = 0, need_w, need_a, W, i, t, side, total, d;
    uint32_t coords[8], v;
    curve_status st;
    curve_range_mode mode = CURVE_RANGE_EXACT;

    need_w = curve_box_intervals_scratch_words(n, bits);
    need_a = curve_box_intervals_scratch_axes(n);
    if (need_w == 0u || need_w > 128u || need_a > 40u) return 0;
    W = curve_key_words_required(n, bits);
    if (W != 1u) return 0;

    if (!oracle_intervals(hilbert, n, bits, bmin, bmax,
                          olo, ohi, ounb, 64, &ocount))
        return 0;

    if (hilbert)
        st = curve_hilbert_box_intervals(n, bits, bmin, bmax,
                                         out_lo, out_hi, out_unb, 64, &count,
                                         mode, scratch, 128, axes, 40);
    else
        st = curve_morton_box_intervals(n, bits, bmin, bmax,
                                        out_lo, out_hi, out_unb, 64, &count,
                                        mode, scratch, 128, axes, 40);
    if (st != CURVE_OK) return 0;
    if (count != ocount) return 0;
    for (i = 0; i < count; i++) {
        if (out_lo[i] != olo[i] || out_unb[i] != ounb[i]) return 0;
        if (!out_unb[i] && out_hi[i] != ohi[i]) return 0;
    }

    /* Count-only. */
    {
        size_t c2 = 99;
        if (hilbert)
            st = curve_hilbert_box_intervals(n, bits, bmin, bmax,
                                             NULL, NULL, NULL, 0, &c2, mode,
                                             scratch, 128, axes, 40);
        else
            st = curve_morton_box_intervals(n, bits, bmin, bmax,
                                            NULL, NULL, NULL, 0, &c2, mode,
                                            scratch, 128, axes, 40);
        if (st != CURVE_OK || c2 != count) return 0;
    }

    /* Every in-box point is covered; no extra points in exact mode. */
    side = 1u << bits;
    total = 1;
    for (d = 0; d < n; d++) total *= side;
    for (t = 0; t < total; t++) {
        uint64_t key;
        int in_box, covered = 0;
        v = t;
        for (d = 0; d < n; d++) {
            coords[d] = v % side;
            v /= side;
        }
        in_box = point_in_box(coords, bmin, bmax, n);
        if (!encode_pt(hilbert, n, bits, coords, &key)) return 0;
        for (i = 0; i < count; i++)
            if (key_in_span(&key, W, &out_lo[i], &out_hi[i], out_unb[i]))
                covered = 1;
        if (in_box && !covered) return 0;
        if (!in_box && covered) return 0;
    }
    return 1;
}

static void test_1d_identity(void)
{
    uint32_t bmin[1] = { 3 }, bmax[1] = { 7 };
    uint64_t scratch[64], lo[4], hi[4];
    int unb[4];
    uint32_t axes[16];
    size_t count = 0, K = 8;
    curve_status st;

    st = curve_morton_box_intervals(1, 8, bmin, bmax, lo, hi, unb, 4, &count,
                                    CURVE_RANGE_EXACT, scratch, 64, axes, 16);
    ASSERT_INT_EQ(st, CURVE_OK, "1D morton");
    ASSERT_INT_EQ(count, 1, "1D one interval");
    ASSERT_INT_EQ(unb[0], 0, "bounded");
    /* 3 and 7 left-aligned in an 8-bit key. */
    ASSERT_INT_EQ(lo[0], (uint64_t)3 << (64 - (int)K), "lower is 3");
    ASSERT_INT_EQ(hi[0], (uint64_t)8 << (64 - (int)K), "exclusive 8");

    st = curve_hilbert_box_intervals(1, 8, bmin, bmax, lo, hi, unb, 4, &count,
                                     CURVE_RANGE_EXACT, scratch, 64, axes, 16);
    ASSERT_INT_EQ(st, CURVE_OK, "1D hilbert");
    ASSERT_INT_EQ(count, 1, "1D hilbert one interval");
    ASSERT_INT_EQ(lo[0], (uint64_t)3 << (64 - (int)K), "hilbert lower 3");
}

static void test_empty_and_full(void)
{
    uint32_t bmin[2] = { 5, 1 }, bmax[2] = { 1, 3 };
    uint64_t scratch[64], lo[4], hi[4];
    int unb[4];
    uint32_t axes[16];
    size_t count = 99;
    curve_status st;

    st = curve_morton_box_intervals(2, 4, bmin, bmax, lo, hi, unb, 4, &count,
                                    CURVE_RANGE_EXACT, scratch, 64, axes, 16);
    ASSERT_INT_EQ(st, CURVE_OK, "inverted box");
    ASSERT_INT_EQ(count, 0, "empty");

    bmin[0] = 0; bmin[1] = 0; bmax[0] = 15; bmax[1] = 15;
    st = curve_morton_box_intervals(2, 4, bmin, bmax, lo, hi, unb, 4, &count,
                                    CURVE_RANGE_EXACT, scratch, 64, axes, 16);
    ASSERT_INT_EQ(st, CURVE_OK, "full domain");
    ASSERT_INT_EQ(count, 1, "full is one interval");
    ASSERT_INT_EQ(unb[0], 1, "full unbounded");
}

static void test_morton_small_grids(void)
{
    unsigned bits;
    int ok = 1;
    for (bits = 1; bits <= 3 && ok; bits++) {
        uint32_t side = 1u << bits, x0, y0, x1, y1;
        for (x0 = 0; x0 < side && ok; x0++)
        for (y0 = 0; y0 < side && ok; y0++)
        for (x1 = x0; x1 < side && ok; x1++)
        for (y1 = y0; y1 < side && ok; y1++) {
            uint32_t bmin[2] = { x0, y0 }, bmax[2] = { x1, y1 };
            if (!check_box(0, 2, bits, bmin, bmax)) ok = 0;
        }
    }
    ASSERT_TRUE(ok, "morton 2D exhaustive bits 1..3");
    {
        uint32_t bmin[3] = { 0, 1, 0 }, bmax[3] = { 2, 3, 1 };
        ASSERT_TRUE(check_box(0, 3, 2, bmin, bmax), "morton 3D bits=2 sample");
    }
}

static void test_hilbert_small_grids(void)
{
    unsigned bits;
    int ok = 1;
    for (bits = 1; bits <= 3 && ok; bits++) {
        uint32_t side = 1u << bits, x0, y0, x1, y1;
        for (x0 = 0; x0 < side && ok; x0++)
        for (y0 = 0; y0 < side && ok; y0++)
        for (x1 = x0; x1 < side && ok; x1++)
        for (y1 = y0; y1 < side && ok; y1++) {
            uint32_t bmin[2] = { x0, y0 }, bmax[2] = { x1, y1 };
            if (!check_box(1, 2, bits, bmin, bmax)) ok = 0;
        }
    }
    ASSERT_TRUE(ok, "hilbert 2D exhaustive bits 1..3");
    {
        uint32_t bmin[3] = { 0, 1, 0 }, bmax[3] = { 2, 3, 1 };
        ASSERT_TRUE(check_box(1, 3, 2, bmin, bmax), "hilbert 3D bits=2 sample");
    }
}

static void test_random_bits4(void)
{
    uint32_t seed = 7, i;
    int ok_m = 1, ok_h = 1;
    for (i = 0; i < 80; i++) {
        uint32_t a = rng_u32(&seed) & 15u, b = rng_u32(&seed) & 15u;
        uint32_t c = rng_u32(&seed) & 15u, d = rng_u32(&seed) & 15u;
        uint32_t bmin[2], bmax[2];
        bmin[0] = a < b ? a : b;
        bmax[0] = a < b ? b : a;
        bmin[1] = c < d ? c : d;
        bmax[1] = c < d ? d : c;
        if (!check_box(0, 2, 4, bmin, bmax)) ok_m = 0;
        if (!check_box(1, 2, 4, bmin, bmax)) ok_h = 0;
    }
    ASSERT_TRUE(ok_m, "morton 2D bits=4 random boxes");
    ASSERT_TRUE(ok_h, "hilbert 2D bits=4 random boxes");
}

static void test_capacity_and_overapprox(void)
{
    uint32_t bmin[2] = { 2, 2 }, bmax[2] = { 13, 13 };
    uint64_t scratch[128], lo[32], hi[32];
    int unb[32];
    uint32_t axes[40];
    size_t exact = 0, clipped = 0, over = 0, t;
    curve_status st;
    uint32_t coords[2];
    int ok;

    st = curve_morton_box_intervals(2, 4, bmin, bmax, lo, hi, unb, 32, &exact,
                                    CURVE_RANGE_EXACT, scratch, 128, axes, 40);
    ASSERT_INT_EQ(st, CURVE_OK, "exact roomy");
    ASSERT_TRUE(exact > 1, "non-trivial box has several intervals");

    st = curve_morton_box_intervals(2, 4, bmin, bmax, lo, hi, unb, 1, &clipped,
                                    CURVE_RANGE_EXACT, scratch, 128, axes, 40);
    ASSERT_INT_EQ(st, CURVE_ENOSPACE, "exact overflow");
    ASSERT_INT_EQ(clipped, exact, "reports full count");

    st = curve_morton_box_intervals(2, 4, bmin, bmax, lo, hi, unb, 1, &over,
                                    CURVE_RANGE_OVERAPPROX, scratch, 128,
                                    axes, 40);
    ASSERT_INT_EQ(st, CURVE_OK, "overapprox");
    ASSERT_TRUE(over <= 1, "overapprox respects capacity");
    ASSERT_INT_EQ(over, 1, "emits one covering span");

    ok = 1;
    for (t = 0; t < 256 && ok; t++) {
        uint64_t key;
        int in_box, covered = 0;
        coords[0] = (uint32_t)(t % 16u);
        coords[1] = (uint32_t)(t / 16u);
        in_box = point_in_box(coords, bmin, bmax, 2);
        encode_pt(0, 2, 4, coords, &key);
        if (key_in_span(&key, 1, &lo[0], &hi[0], unb[0])) covered = 1;
        if (in_box && !covered) ok = 0;
    }
    ASSERT_TRUE(ok, "overapprox has no false negatives");
}

static void test_range_failures(void)
{
    uint32_t bmin[2] = { 0, 0 }, bmax[2] = { 1, 1 };
    uint64_t scratch[64], lo[2], hi[2];
    int unb[2];
    uint32_t axes[16];
    size_t count = 0;

    ASSERT_INT_EQ(curve_morton_box_intervals(9, 4, bmin, bmax, lo, hi, unb, 2,
                                             &count, CURVE_RANGE_EXACT,
                                             scratch, 64, axes, 16),
                  CURVE_EUNSUPPORTED, "n=9");
    ASSERT_INT_EQ(curve_morton_box_intervals(2, 4, NULL, bmax, lo, hi, unb, 2,
                                             &count, CURVE_RANGE_EXACT,
                                             scratch, 64, axes, 16),
                  CURVE_EINVAL, "null box");
    ASSERT_INT_EQ(curve_morton_box_intervals(2, 4, bmin, bmax, lo, hi, unb, 2,
                                             NULL, CURVE_RANGE_EXACT,
                                             scratch, 64, axes, 16),
                  CURVE_EINVAL, "null count");
    ASSERT_INT_EQ(curve_hilbert_box_intervals(2, 4, bmin, bmax, lo, hi, unb, 2,
                                              &count, CURVE_RANGE_EXACT,
                                              scratch, 4, axes, 16),
                  CURVE_ENOSPACE, "short scratch");
}

static void test_prefix_dir(void)
{
    enum { N = 16 };
    uint64_t keys[N * N], scratch[4 * (N * N + 1)], pr_scratch[4];
    uint32_t coords[2], sc[2], mins[512], maxs[512], bmin[2], bmax[2];
    curve_prefix_node nodes[512];
    size_t n_nodes = 0, i, x, y, begin = 0, end = 0, n_ranges;
    size_t rb[32], re[32];
    curve_status st;
    int ok;

    i = 0;
    for (y = 0; y < N; y++) {
        for (x = 0; x < N; x++) {
            coords[0] = (uint32_t)x;
            coords[1] = (uint32_t)y;
            curve_hilbertnd_encode_words(coords, 2, 4, &keys[i], 1, sc, 2);
            i++;
        }
    }
    g_sort_words = 1;
    qsort(keys, (size_t)N * N, sizeof(uint64_t), cmp_keys);

    st = curve_prefix_dir_build(keys, (size_t)N * N, 1, 8, 4,
                                nodes, 512, &n_nodes, scratch,
                                4 * (N * N + 1));
    ASSERT_INT_EQ(st, CURVE_OK, "build dir");
    ASSERT_TRUE(n_nodes >= 1 && n_nodes <= 512, "node count");

    st = curve_prefix_dir_fill_bounds(nodes, n_nodes, keys, 1, 2, 4, 1,
                                      mins, maxs, sc, 2);
    ASSERT_INT_EQ(st, CURVE_OK, "fill bounds");

    ok = 1;
    for (i = 0; i < (size_t)N * N && ok; i++) {
        size_t b2, e2;
        st = curve_prefix_dir_prefix(nodes, n_nodes, keys, 1, 8,
                                     &keys[i], 4, pr_scratch, 4, &begin, &end);
        if (st != CURVE_OK) ok = 0;
        curve_key_prefix_range(keys, (size_t)N * N, 1, &keys[i], 8, 4,
                               pr_scratch, 4, &b2, &e2);
        if (begin != b2 || end != e2) ok = 0;
    }
    ASSERT_TRUE(ok, "dir prefix lookup matches prefix_range");

    bmin[0] = 3; bmin[1] = 3; bmax[0] = 10; bmax[1] = 12;
    n_ranges = 0;
    st = curve_prefix_dir_box(nodes, n_nodes, mins, maxs, 2, bmin, bmax,
                              keys, 1, 4, 1, sc, 2, scratch, n_nodes,
                              rb, re, 32, &n_ranges);
    ASSERT_INT_EQ(st, CURVE_OK, "dir box");
    ASSERT_TRUE(n_ranges >= 1, "at least one index range");

    ok = 1;
    for (i = 0; i < (size_t)N * N && ok; i++) {
        uint32_t dec[2];
        int in_box, in_range = 0;
        size_t r;
        curve_hilbertnd_decode_words(&keys[i], 1, 2, 4, dec, 2, sc, 2);
        in_box = point_in_box(dec, bmin, bmax, 2);
        for (r = 0; r < n_ranges; r++)
            if (i >= rb[r] && i < re[r]) in_range = 1;
        if (in_box && !in_range) ok = 0;
        if (!in_box && in_range) ok = 0;
    }
    ASSERT_TRUE(ok, "exact dir box matches brute");

    {
        uint64_t unsorted[3] = { 2, 0, 1 };
        size_t nn = 0;
        ASSERT_INT_EQ(curve_prefix_dir_build(unsorted, 3, 1, 8, 1, nodes, 16,
                                             &nn, scratch, 16),
                      CURVE_EINVAL, "unsorted rejected");
    }
}

int main(void)
{
    RUN_TEST(test_1d_identity);
    RUN_TEST(test_empty_and_full);
    RUN_TEST(test_morton_small_grids);
    RUN_TEST(test_hilbert_small_grids);
    RUN_TEST(test_random_bits4);
    RUN_TEST(test_capacity_and_overapprox);
    RUN_TEST(test_range_failures);
    RUN_TEST(test_prefix_dir);
    TEST_SUMMARY();
}
