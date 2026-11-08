#include "curves.h"
#include "test_common.h"
#include "ref_curves.h"
#include <stdlib.h>

/* -------------------------------------------------------------------------- */
/* Gray codes                                                                 */
/* -------------------------------------------------------------------------- */

static void test_gray_inverse_and_adjacency(void)
{
    int ok = 1;
    for (uint64_t i = 0; i < 1000000; ++i) {
        if (curve_gray_decode(curve_gray_encode(i)) != i) { ok = 0; break; }
    }
    ASSERT_TRUE(ok, "gray encode/decode inverse over 1M values");

    ok = 1;
    for (uint64_t i = 0; i < 999999; ++i) {
        uint64_t diff = curve_gray_encode(i) ^ curve_gray_encode(i + 1);
        if (diff == 0 || (diff & (diff - 1)) != 0) { ok = 0; break; }
    }
    ASSERT_TRUE(ok, "adjacent Gray codes differ by one bit");

    ASSERT_INT_EQ(curve_gray_encode(0), 0, "gray(0) == 0");
    ASSERT_INT_EQ(curve_gray_encode(1), 1, "gray(1) == 1");
    ASSERT_INT_EQ(curve_gray_encode(2), 3, "gray(2) == 3");
}

/* -------------------------------------------------------------------------- */
/* Morton / Z-order                                                           */
/* -------------------------------------------------------------------------- */

static void test_morton2_roundtrip(void)
{
    int ok = 1;
    for (uint32_t x = 0; x < 256; ++x)
        for (uint32_t y = 0; y < 256; ++y) {
            uint32_t d = curve_morton2(x, y);
            uint32_t rx, ry;
            curve_morton2_decode(d, &rx, &ry);
            if (rx != x || ry != y) { ok = 0; break; }
        }
    ASSERT_TRUE(ok, "morton2 encode/decode roundtrip");

    /* known answers: x on even bits, y on odd bits */
    ASSERT_INT_EQ(curve_morton2(1, 0), 1, "morton2(1,0) == 1");
    ASSERT_INT_EQ(curve_morton2(0, 1), 2, "morton2(0,1) == 2");
    ASSERT_INT_EQ(curve_morton2(3, 2), 13, "morton2(3,2) == 0b1101");
    ASSERT_INT_EQ(curve_morton2(0xFFFF, 0xFFFF), 0xFFFFFFFFu, "morton2 max");

    ok = 1;
    for (uint32_t x = 0; x < 256 && ok; ++x)
        for (uint32_t y = 0; y < 256 && ok; ++y) {
            uint32_t d = curve_morton2(x, y);
            uint32_t rx, ry;
            if (d != ref_morton2(x, y)) { ok = 0; break; }
            ref_morton2_decode(d, &rx, &ry);
            if (rx != x || ry != y) { ok = 0; break; }
        }
    ASSERT_TRUE(ok, "morton2 matches nested-loop reference");
}

static void test_morton3_roundtrip(void)
{
    int ok = 1;
    for (uint32_t x = 0; x < 32; ++x)
        for (uint32_t y = 0; y < 32; ++y)
            for (uint32_t z = 0; z < 32; ++z) {
                uint64_t d = curve_morton3(x, y, z);
                uint32_t rx, ry, rz;
                curve_morton3_decode(d, &rx, &ry, &rz);
                if (rx != x || ry != y || rz != z) { ok = 0; goto done; }
            }
done:
    ASSERT_TRUE(ok, "morton3 encode/decode roundtrip");

    ASSERT_INT_EQ(curve_morton3(1, 0, 0), 1, "morton3(1,0,0) == 1");
    ASSERT_INT_EQ(curve_morton3(0, 1, 0), 2, "morton3(0,1,0) == 2");
    ASSERT_INT_EQ(curve_morton3(0, 0, 1), 4, "morton3(0,0,1) == 4");

    ok = 1;
    for (uint32_t x = 0; x < 16 && ok; ++x)
        for (uint32_t y = 0; y < 16 && ok; ++y)
            for (uint32_t z = 0; z < 16 && ok; ++z) {
                uint64_t d = curve_morton3(x, y, z);
                uint32_t rx, ry, rz;
                if (d != ref_morton3(x, y, z)) { ok = 0; break; }
                ref_morton3_decode(d, &rx, &ry, &rz);
                if (rx != x || ry != y || rz != z) { ok = 0; break; }
            }
    ASSERT_TRUE(ok, "morton3 matches nested-loop reference, 16^3");
}

static void test_morton2_64(void)
{
    uint64_t x, y, d;
    uint32_t rx, ry;
    int ok = 1;

    /* consistency with 32-bit version on the low range */
    for (x = 0; x < 512 && ok; ++x)
        for (y = 0; y < 512 && ok; ++y) {
            uint64_t got = curve_morton2_64((uint32_t)x, (uint32_t)y);
            if (got != curve_morton2((uint32_t)x, (uint32_t)y) ||
                got != ref_morton2_64((uint32_t)x, (uint32_t)y))
                ok = 0;
        }
    ASSERT_TRUE(ok, "morton2_64 consistent with morton2 on low range");

    /* roundtrip including high bits */
    ok = 1;
    for (int i = 0; i < 100000 && ok; ++i) {
        x = (uint64_t)i * 2654435761u;
        y = (uint64_t)i * 2246822519u + 12345u;
        x &= 0xFFFFFFFFull; y &= 0xFFFFFFFFull;
        d = curve_morton2_64((uint32_t)x, (uint32_t)y);
        curve_morton2_64_decode(d, &rx, &ry);
        if (rx != (uint32_t)x || ry != (uint32_t)y) ok = 0;
    }
    ASSERT_TRUE(ok, "morton2_64 roundtrip with full 32-bit axes");

    ASSERT_INT_EQ(curve_morton2_64(0xFFFFFFFFu, 0xFFFFFFFFu),
                  0xFFFFFFFFFFFFFFFFull, "morton2_64 max");
    ASSERT_INT_EQ(curve_morton2_64(0x80000000u, 0), 0x4000000000000000ull,
                  "morton2_64 high bit lands at bit 62");
}

/* -------------------------------------------------------------------------- */
/* Hilbert 2D (32-bit and 64-bit)                                             */
/* -------------------------------------------------------------------------- */

static void test_hilbert2_known_vectors(void)
{
    /* Burkardt hilbert_curve reference table (8x8 grid, bits=3) */
    uint32_t x, y;
    curve_hilbert_d2xy(0, 3, &x, &y);
    ASSERT_TRUE(x == 0 && y == 0, "hilbert2 d=0 -> (0,0)");
    curve_hilbert_d2xy(1, 3, &x, &y);
    ASSERT_TRUE(x == 0 && y == 1, "hilbert2 d=1 -> (0,1)");
    curve_hilbert_d2xy(2, 3, &x, &y);
    ASSERT_TRUE(x == 1 && y == 1, "hilbert2 d=2 -> (1,1)");
    curve_hilbert_d2xy(3, 3, &x, &y);
    ASSERT_TRUE(x == 1 && y == 0, "hilbert2 d=3 -> (1,0)");
    curve_hilbert_d2xy(63, 3, &x, &y);
    ASSERT_TRUE(x == 7 && y == 0, "hilbert2 d=63 -> (7,0)");

    ASSERT_INT_EQ(curve_hilbert_xy2d(0, 7, 3), 21, "hilbert2 (0,7) -> 21");
    ASSERT_INT_EQ(curve_hilbert_xy2d(7, 0, 3), 63, "hilbert2 (7,0) -> 63");
    ASSERT_INT_EQ(curve_hilbert_xy2d(3, 3, 3), 10, "hilbert2 (3,3) -> 10");
    ASSERT_INT_EQ(curve_hilbert_xy2d(4, 4, 3), 32, "hilbert2 (4,4) -> 32");
}

static void test_hilbert2_bijective_all_bits(void)
{
    int ok = 1;
    for (int bits = 1; bits <= 8 && ok; ++bits) {
        int N = 1 << bits;
        int total = N * N;
        uint8_t *seen = (uint8_t*)calloc((size_t)total, 1);
        if (!seen) { ASSERT_TRUE(0, "malloc"); return; }
        for (int x = 0; x < N && ok; ++x)
            for (int y = 0; y < N && ok; ++y) {
                uint32_t d = curve_hilbert_xy2d((uint32_t)x, (uint32_t)y, bits);
                uint32_t rx, ry;
                if (d >= (uint32_t)total || seen[d]) { ok = 0; break; }
                seen[d] = 1;
                curve_hilbert_d2xy(d, bits, &rx, &ry);
                if (rx != (uint32_t)x || ry != (uint32_t)y) { ok = 0; break; }
                /* cross-check against slow reference */
                if (ref_hilbert_xy2d((uint64_t)x, (uint64_t)y, bits) != d) {
                    ok = 0; break;
                }
                {
                    uint64_t ux, uy;
                    ref_hilbert_d2xy(d, bits, &ux, &uy);
                    if (ux != (uint64_t)x || uy != (uint64_t)y) {
                        ok = 0; break;
                    }
                }
            }
        free(seen);
    }
    ASSERT_TRUE(ok, "hilbert2 bijective + inverse + matches reference, bits 1..8");
}

static void test_hilbert2_adjacency(void)
{
    int ok = 1;
    for (int bits = 1; bits <= 5 && ok; ++bits) {
        uint32_t N = 1u << bits;
        uint64_t total = (uint64_t)N * N;
        uint32_t px = 0, py = 0;
        for (uint64_t d = 0; d < total; ++d) {
            uint32_t x, y;
            curve_hilbert_d2xy((uint32_t)d, bits, &x, &y);
            if (d > 0) {
                uint32_t dx = x > px ? x - px : px - x;
                uint32_t dy = y > py ? y - py : py - y;
                if (dx + dy != 1) { ok = 0; break; }
            }
            px = x;
            py = y;
        }
    }
    ASSERT_TRUE(ok, "hilbert2 consecutive cells are edge-adjacent, bits 1..5");
}

static void test_hilbert2_sampled_16bit(void)
{
    int ok = 1;
    for (uint32_t i = 0; i < 200000 && ok; ++i) {
        uint32_t x = (i * 2654435761u) & 0xFFFFu;
        uint32_t y = (i * 2246822519u + 77u) & 0xFFFFu;
        uint32_t d = curve_hilbert_xy2d(x, y, 16);
        uint32_t rx, ry;
        curve_hilbert_d2xy(d, 16, &rx, &ry);
        if (rx != x || ry != y) ok = 0;
    }
    ASSERT_TRUE(ok, "hilbert2 16-bit sampled roundtrip");
}

static void test_hilbert2_64(void)
{
    int ok = 1;

    /* consistency with 32-bit version, exhaustive at low bits */
    for (int bits = 1; bits <= 8 && ok; ++bits) {
        uint32_t N = 1u << bits;
        for (uint32_t x = 0; x < N && ok; ++x)
            for (uint32_t y = 0; y < N && ok; ++y) {
                uint64_t d64 = curve_hilbert2_64_xy2d(x, y, bits);
                uint32_t d32 = curve_hilbert_xy2d(x, y, bits);
                uint32_t rx, ry;
                if (d64 != d32) { ok = 0; break; }
                curve_hilbert2_64_d2xy(d64, bits, &rx, &ry);
                if (rx != x || ry != y) { ok = 0; break; }
                if (ref_hilbert_xy2d(x, y, bits) != d64) { ok = 0; break; }
            }
    }
    ASSERT_TRUE(ok, "hilbert2_64 matches 32-bit version and reference, bits 1..8");

    /* sampled roundtrip at high bit counts */
    ok = 1;
    for (int bits = 17; bits <= 32 && ok; ++bits) {
        for (uint32_t i = 0; i < 20000 && ok; ++i) {
            uint32_t x = (i * 2654435761u);
            uint32_t y = (i * 2246822519u + 99u);
            if (bits < 32) { x >>= (32 - bits); y >>= (32 - bits); }
            uint64_t d = curve_hilbert2_64_xy2d(x, y, bits);
            uint32_t rx, ry;
            curve_hilbert2_64_d2xy(d, bits, &rx, &ry);
            if (rx != x || ry != y) ok = 0;
        }
    }
    ASSERT_TRUE(ok, "hilbert2_64 sampled roundtrip, bits 17..32");
}

/* -------------------------------------------------------------------------- */
/* Hilbert 3D                                                                 */
/* -------------------------------------------------------------------------- */

static int hilbert3_check_bits(int bits)
{
    int N = 1 << bits;
    int total = N * N * N;
    int ok = 1;
    uint8_t *seen = (uint8_t*)calloc((size_t)total, 1);
    if (!seen) return 0;
    for (int x = 0; x < N && ok; ++x)
        for (int y = 0; y < N && ok; ++y)
            for (int z = 0; z < N && ok; ++z) {
                uint32_t d = curve_hilbert3_xyz2d((uint32_t)x, (uint32_t)y,
                                                  (uint32_t)z, bits);
                uint32_t rx, ry, rz;
                if (d >= (uint32_t)total || seen[d]) { ok = 0; break; }
                seen[d] = 1;
                curve_hilbert3_d2xyz(d, bits, &rx, &ry, &rz);
                if (rx != (uint32_t)x || ry != (uint32_t)y || rz != (uint32_t)z) {
                    ok = 0; break;
                }
            }
    free(seen);
    return ok;
}

static void test_hilbert3_bijective(void)
{
    ASSERT_TRUE(hilbert3_check_bits(1), "hilbert3 bijective at 1 bit");
    ASSERT_TRUE(hilbert3_check_bits(4), "hilbert3 bijective at 4 bits");
    ASSERT_TRUE(hilbert3_check_bits(5), "hilbert3 bijective at 5 bits");

    /* sampled at 10 bits */
    int ok = 1;
    for (uint32_t i = 0; i < 100000 && ok; ++i) {
        uint32_t x = (i * 2654435761u) & 1023u;
        uint32_t y = (i * 2246822519u + 7u) & 1023u;
        uint32_t z = (i * 3266489917u + 13u) & 1023u;
        uint32_t d = curve_hilbert3_xyz2d(x, y, z, 10);
        uint32_t rx, ry, rz;
        curve_hilbert3_d2xyz(d, 10, &rx, &ry, &rz);
        if (rx != x || ry != y || rz != z) ok = 0;
    }
    ASSERT_TRUE(ok, "hilbert3 10-bit sampled roundtrip");
}

/* -------------------------------------------------------------------------- */
/* Hilbert n-D                                                                */
/* -------------------------------------------------------------------------- */

static void test_hilbertnd_roundtrip(void)
{
    int ok = 1;
    for (int n = 2; n <= 8 && ok; ++n) {
        for (int bits = 1; bits <= 8 && ok; ++bits) {
            if (n * bits > 24) break;   /* keep sample counts sane */
            uint64_t total = 1ull << (n * bits);
            uint64_t step = total > 4096 ? total / 4096 : 1;
            for (uint64_t d = 0; d < total && ok; d += step) {
                uint32_t coords[8], back[8];
                for (int i = 0; i < n; ++i) coords[i] = 0;
                curve_hilbertnd_d2xy(d, coords, n, bits);
                for (int i = 0; i < n; ++i)
                    if (coords[i] >= (1u << bits)) { ok = 0; break; }
                if (!ok) break;
                if (curve_hilbertnd_xy2d(coords, n, bits) != d) { ok = 0; break; }
                /* mutual-inverse sanity via double application */
                curve_hilbertnd_d2xy(curve_hilbertnd_xy2d(coords, n, bits),
                                     back, n, bits);
                for (int i = 0; i < n; ++i)
                    if (back[i] != coords[i]) { ok = 0; break; }
            }
        }
    }
    ASSERT_TRUE(ok, "hilbertnd roundtrip for n=2..8, bits=1..8");
}

static void test_hilbertnd_matches_2d(void)
{
    int ok = 1;
    for (int bits = 1; bits <= 8 && ok; ++bits) {
        uint32_t N = 1u << bits;
        for (uint32_t x = 0; x < N && ok; ++x)
            for (uint32_t y = 0; y < N && ok; ++y) {
                uint32_t coords[2] = { x, y };
                if (curve_hilbertnd_xy2d(coords, 2, bits) !=
                    curve_hilbert_xy2d(x, y, bits))
                    ok = 0;
            }
    }
    ASSERT_TRUE(ok, "hilbertnd n=2 matches curve_hilbert_xy2d");
}

/* One Hilbert level is n_dims bits of the index. Points sharing that prefix
 * occupy an axis-aligned 2^(bits-1) hypercube. */
static void test_hilbertnd_prefix_cells(void)
{
    int ok = 1;
    int bits = 4;
    uint32_t N = 1u << bits;
    uint32_t half = N / 2u;
    uint64_t total = (uint64_t)N * N;
    uint64_t cell = total / 4u; /* 2-bit prefix → 4 cells */
    uint64_t c;
    for (c = 0; c < 4 && ok; ++c) {
        uint32_t minx = 0xFFFFFFFFu, miny = 0xFFFFFFFFu, maxx = 0, maxy = 0;
        uint64_t d;
        for (d = c * cell; d < (c + 1) * cell; ++d) {
            uint32_t coords[2];
            curve_hilbertnd_d2xy(d, coords, 2, bits);
            if (coords[0] < minx) minx = coords[0];
            if (coords[1] < miny) miny = coords[1];
            if (coords[0] > maxx) maxx = coords[0];
            if (coords[1] > maxy) maxy = coords[1];
        }
        if (maxx - minx + 1u != half || maxy - miny + 1u != half)
            ok = 0;
    }
    ASSERT_TRUE(ok, "hilbertnd n=2 bits=4: 2-bit prefixes are 8x8 squares");
}

static void test_ref_pack_roundtrip(void)
{
    uint32_t X[4] = { 1u, 2u, 3u, 4u }, Y[4];
    uint64_t h, words[2];
    int i, ok = 1;
    h = ref_pack_axes_u64(X, 4, 4);
    ref_unpack_axes_u64(h, 4, 4, Y);
    for (i = 0; i < 4; i++)
        if (Y[i] != X[i]) ok = 0;
    ref_pack_axes_words(X, 4, 4, words, 2);
    ASSERT_TRUE(ok && words[0] == (h << (64 - 16)) && words[1] == 0,
                "ref pack u64 matches left-aligned words");
    ref_unpack_axes_words(words, 4, 4, Y);
    ok = 1;
    for (i = 0; i < 4; i++)
        if (Y[i] != X[i]) ok = 0;
    ASSERT_TRUE(ok, "ref pack words roundtrip");
}

/* -------------------------------------------------------------------------- */
/* Moore curve                                                                */
/* -------------------------------------------------------------------------- */

static void test_moore_known_path(void)
{
    /* order-2 path derived from the Moore L-system */
    static const uint8_t path[16][2] = {
        {1,0}, {0,0}, {0,1}, {1,1}, {1,2}, {0,2}, {0,3}, {1,3},
        {2,3}, {3,3}, {3,2}, {2,2}, {2,1}, {3,1}, {3,0}, {2,0}
    };
    int ok = 1;
    for (uint32_t d = 0; d < 16 && ok; ++d) {
        uint32_t x, y;
        curve_moore_d2xy(d, 2, &x, &y);
        if (x != path[d][0] || y != path[d][1]) ok = 0;
        if (curve_moore_xy2d(path[d][0], path[d][1], 2) != d) ok = 0;
    }
    ASSERT_TRUE(ok, "moore order-2 matches L-system-derived path");
}

static void test_moore_bijective_continuous(void)
{
    int ok = 1, closed = 1, continuous = 1;
    for (int bits = 1; bits <= 8 && ok; ++bits) {
        uint32_t N = 1u << bits;
        uint64_t total = (uint64_t)N * N;
        uint8_t *seen = (uint8_t*)calloc((size_t)total, 1);
        uint32_t fx = 0, fy = 0, px = 0, py = 0;
        if (!seen) { ASSERT_TRUE(0, "malloc"); return; }
        for (uint64_t d = 0; d < total; ++d) {
            uint32_t x, y;
            curve_moore_d2xy((uint32_t)d, bits, &x, &y);
            if (x >= N || y >= N || seen[(uint64_t)y * N + x]) { ok = 0; break; }
            seen[(uint64_t)y * N + x] = 1;
            if (curve_moore_xy2d(x, y, bits) != d) { ok = 0; break; }
            if (d == 0) { fx = x; fy = y; }
            else {
                uint32_t dx = x > px ? x - px : px - x;
                uint32_t dy = y > py ? y - py : py - y;
                if (dx + dy != 1) continuous = 0;
            }
            px = x; py = y;
        }
        {
            uint32_t dx = fx > px ? fx - px : px - fx;
            uint32_t dy = fy > py ? fy - py : py - fy;
            if (dx + dy != 1) closed = 0;
        }
        free(seen);
    }
    ASSERT_TRUE(ok, "moore bijective + inverse, bits 1..8");
    ASSERT_TRUE(continuous, "moore consecutive points edge-adjacent");
    ASSERT_TRUE(closed, "moore curve closes into a loop");
}

/* -------------------------------------------------------------------------- */
/* Sierpiński ordering key                                                    */
/* -------------------------------------------------------------------------- */

static int cmp_u64(const void *a, const void *b)
{
    uint64_t va = *(const uint64_t*)a, vb = *(const uint64_t*)b;
    return (va > vb) - (va < vb);
}

typedef struct { uint64_t k; uint32_t x, y; } sierp_ent;

static int cmp_sierp_ent(const void *a, const void *b)
{
    uint64_t va = ((const sierp_ent*)a)->k, vb = ((const sierp_ent*)b)->k;
    return (va > vb) - (va < vb);
}

static void test_sierpinski_ordering(void)
{
    int injective = 1, continuous = 1, closed = 1;

    ASSERT_INT_EQ(curve_sierpinski_xy2d(0, 0, 4), 0, "sierpinski(0,0) == 0");

    for (int bits = 1; bits <= 7 && injective; ++bits) {
        uint32_t N = 1u << bits;
        uint64_t total = (uint64_t)N * N;
        /* key space is 2*total entries (non-dense) */
        uint64_t *keys = (uint64_t*)malloc((size_t)total * sizeof(uint64_t));
        if (!keys) { ASSERT_TRUE(0, "malloc"); return; }
        for (uint32_t y = 0; y < N; ++y)
            for (uint32_t x = 0; x < N; ++x) {
                uint64_t k = curve_sierpinski_xy2d(x, y, bits);
                if (k >= 2 * total) { injective = 0; }
                keys[(uint64_t)y * N + x] = k;
            }
        /* injectivity via sort */
        qsort(keys, (size_t)total, sizeof(uint64_t), cmp_u64);
        for (uint64_t i = 1; i < total; ++i)
            if (keys[i] == keys[i - 1]) { injective = 0; break; }
        free(keys);

        /* 8-neighbor continuity of the induced order */
        if (injective && bits <= 5) {
            uint64_t *order = (uint64_t*)malloc((size_t)total * sizeof(uint64_t));
            uint32_t *cx = (uint32_t*)malloc((size_t)total * sizeof(uint32_t));
            uint32_t *cy = (uint32_t*)malloc((size_t)total * sizeof(uint32_t));
            if (!order || !cx || !cy) { ASSERT_TRUE(0, "malloc"); return; }
            for (uint32_t y = 0; y < N; ++y)
                for (uint32_t x = 0; x < N; ++x) {
                    uint64_t i = (uint64_t)y * N + x;
                    order[i] = curve_sierpinski_xy2d(x, y, bits);
                    cx[i] = x; cy[i] = y;
                }
            sierp_ent *ents = (sierp_ent*)malloc((size_t)total * sizeof(sierp_ent));
            if (!ents) { ASSERT_TRUE(0, "malloc"); return; }
            for (uint64_t i = 0; i < total; ++i) {
                ents[i].k = order[i]; ents[i].x = cx[i]; ents[i].y = cy[i];
            }
            qsort(ents, (size_t)total, sizeof(sierp_ent), cmp_sierp_ent);
            for (uint64_t i = 1; i < total; ++i) {
                uint32_t dx = ents[i].x > ents[i-1].x ? ents[i].x - ents[i-1].x
                                                      : ents[i-1].x - ents[i].x;
                uint32_t dy = ents[i].y > ents[i-1].y ? ents[i].y - ents[i-1].y
                                                      : ents[i-1].y - ents[i].y;
                if (dx > 1 || dy > 1 || (dx == 0 && dy == 0)) continuous = 0;
            }
            {
                uint32_t dx = ents[0].x > ents[total-1].x ? ents[0].x - ents[total-1].x
                                                          : ents[total-1].x - ents[0].x;
                uint32_t dy = ents[0].y > ents[total-1].y ? ents[0].y - ents[total-1].y
                                                          : ents[total-1].y - ents[0].y;
                if (dx > 1 || dy > 1) closed = 0;
            }
            free(ents); free(order); free(cx); free(cy);
        }
    }
    ASSERT_TRUE(injective, "sierpinski key injective, in range");
    ASSERT_TRUE(continuous, "sierpinski order 8-neighbor continuous");
    ASSERT_TRUE(closed, "sierpinski order closes into a loop");
}

/* -------------------------------------------------------------------------- */
/* Windowed helpers                                                           */
/* -------------------------------------------------------------------------- */

static void test_windowed_matches_discrete(void)
{
    const double w[4] = { -2.0, 3.0, 7.5, 11.25 };
    int ok = 1;
    for (int i = 0; i <= 40 && ok; ++i)
        for (int j = 0; j <= 40 && ok; ++j) {
            double x = w[0] + (w[2] - w[0]) * (i / 40.0);
            double y = w[1] + (w[3] - w[1]) * (j / 40.0);
            uint32_t xi = (uint32_t)((x - w[0]) / (w[2] - w[0]) * 0xFFFF);
            uint32_t yi = (uint32_t)((y - w[1]) / (w[3] - w[1]) * 0xFFFF);
            if (curve_hilbert_window(x, y, w) != curve_hilbert_xy2d(xi, yi, 16))
                ok = 0;
            if (curve_z_window(x, y, w) != curve_morton2(xi, yi))
                ok = 0;
        }
    ASSERT_TRUE(ok, "windowed helpers match discrete grid over window");

    /* degenerate window (both axes zero-size) maps to 0 */
    const double wdeg[4] = { 1.0, 1.0, 1.0, 1.0 };
    ASSERT_INT_EQ(curve_z_window(1.0, 2.0, wdeg), 0, "degenerate window");
    ASSERT_INT_EQ(curve_hilbert_window(1.0, 2.0, wdeg),
                  curve_hilbert_xy2d(0, 0, 16), "degenerate window hilbert");
}

/* -------------------------------------------------------------------------- */
/* Gilbert curves                                                             */
/* -------------------------------------------------------------------------- */

typedef struct {
    int w, h, d;
    int *seen;
    long count;
    int oob, discontinuous;
    int px, py, pz;
} gilbert_ctx;

static void gilbert2_mark(int x, int y, void *user)
{
    gilbert_ctx *c = (gilbert_ctx*)user;
    if (x < 0 || y < 0 || x >= c->w || y >= c->h) { c->oob = 1; return; }
    c->seen[y * c->w + x]++;
    if (c->count > 0) {
        int dx = x > c->px ? x - c->px : c->px - x;
        int dy = y > c->py ? y - c->py : c->py - y;
        if (dx + dy != 1) c->discontinuous = 1;
    }
    c->px = x; c->py = y;
    c->count++;
}

static void gilbert3_mark(int x, int y, int z, void *user)
{
    gilbert_ctx *c = (gilbert_ctx*)user;
    if (x < 0 || y < 0 || z < 0 || x >= c->w || y >= c->h || z >= c->d) {
        c->oob = 1;
        return;
    }
    c->seen[(z * c->h + y) * c->w + x]++;
    if (c->count > 0) {
        int dx = x > c->px ? x - c->px : c->px - x;
        int dy = y > c->py ? y - c->py : c->py - y;
        int dz = z > c->pz ? z - c->pz : c->pz - z;
        if (dx + dy + dz != 1) c->discontinuous = 1;
    }
    c->px = x; c->py = y; c->pz = z;
    c->count++;
}

static int gilbert2_check(int w, int h)
{
    int ok;
    int *seen = (int*)calloc((size_t)(w * h), sizeof(int));
    gilbert_ctx ctx = { w, h, 0, seen, 0, 0, 0, 0, 0, 0 };
    if (!seen) return 0;
    curve_gilbert2(w, h, gilbert2_mark, &ctx);
    ok = (ctx.count == (long)w * h) && !ctx.oob && !ctx.discontinuous;
    if (ok)
        for (int i = 0; i < w * h; ++i)
            if (seen[i] != 1) { ok = 0; break; }
    free(seen);
    return ok;
}

static int gilbert3_check(int w, int h, int d)
{
    int ok;
    /* upstream guarantees face-adjacency only for even dimensions */
    int need_continuity = !(w & 1) && !(h & 1) && !(d & 1);
    long n = (long)w * h * d;
    int *seen = (int*)calloc((size_t)n, sizeof(int));
    gilbert_ctx ctx = { w, h, d, seen, 0, 0, 0, 0, 0, 0 };
    if (!seen) return 0;
    curve_gilbert3(w, h, d, gilbert3_mark, &ctx);
    ok = (ctx.count == n) && !ctx.oob && !(need_continuity && ctx.discontinuous);
    if (ok)
        for (long i = 0; i < n; ++i)
            if (seen[i] != 1) { ok = 0; break; }
    free(seen);
    return ok;
}

static void test_gilbert2_exact_cover_continuous(void)
{
    static const int shapes[][2] = {
        {1, 1}, {1, 7}, {7, 1}, {2, 2}, {5, 7}, {8, 8}, {13, 7}, {16, 16},
        {2, 3}, {3, 2}, {100, 63}
    };
    int ok = 1;
    for (size_t i = 0; i < sizeof(shapes) / sizeof(shapes[0]); ++i)
        if (!gilbert2_check(shapes[i][0], shapes[i][1])) {
            ok = 0;
            fprintf(stderr, "    gilbert2 failed for %dx%d\n",
                    shapes[i][0], shapes[i][1]);
        }
    ASSERT_TRUE(ok, "gilbert2 exact cover + continuity, incl. degenerate shapes");
}

static void test_gilbert3_exact_cover_continuous(void)
{
    static const int shapes[][3] = {
        {1, 1, 1}, {1, 1, 9}, {1, 4, 3}, {4, 3, 5}, {2, 2, 2}, {5, 3, 3},
        {1, 5, 7}, {8, 8, 4}, {6, 4, 4}, {4, 4, 4}
    };
    int ok = 1;
    for (size_t i = 0; i < sizeof(shapes) / sizeof(shapes[0]); ++i)
        if (!gilbert3_check(shapes[i][0], shapes[i][1], shapes[i][2])) {
            ok = 0;
            fprintf(stderr, "    gilbert3 failed for %dx%dx%d\n",
                    shapes[i][0], shapes[i][1], shapes[i][2]);
        }
    ASSERT_TRUE(ok, "gilbert3 exact cover + continuity, incl. degenerate boxes");
}

/* -------------------------------------------------------------------------- */
/* Total-function clamping                                                    */
/* -------------------------------------------------------------------------- */

static void test_clamping(void)
{
    uint32_t x, y;
    /* bits out of range must clamp, never crash */
    ASSERT_INT_EQ(curve_hilbert_xy2d(1, 1, 0), curve_hilbert_xy2d(1, 1, 1),
                  "hilbert2 bits=0 clamps to 1");
    ASSERT_INT_EQ(curve_hilbert_xy2d(5, 5, 99), curve_hilbert_xy2d(5, 5, 16),
                  "hilbert2 bits=99 clamps to 16");
    ASSERT_INT_EQ(curve_hilbert3_xyz2d(1, 1, 1, -3),
                  curve_hilbert3_xyz2d(1, 1, 1, 1), "hilbert3 bits<1 clamps");
    ASSERT_INT_EQ(curve_hilbert2_64_xy2d(5, 5, 0),
                  curve_hilbert2_64_xy2d(5, 5, 1), "hilbert2_64 bits=0 clamps");
    ASSERT_INT_EQ(curve_hilbert2_64_xy2d(5, 5, 77),
                  curve_hilbert2_64_xy2d(5, 5, 32), "hilbert2_64 bits=77 clamps");
    ASSERT_INT_EQ(curve_moore_xy2d(3, 3, 0), curve_moore_xy2d(3, 3, 1),
                  "moore bits=0 clamps");
    ASSERT_INT_EQ(curve_moore_xy2d(3, 3, 50), curve_moore_xy2d(3, 3, 16),
                  "moore bits=50 clamps");
    ASSERT_INT_EQ(curve_sierpinski_xy2d(3, 3, 0),
                  curve_sierpinski_xy2d(3, 3, 1), "sierpinski bits=0 clamps");

    {
        /* bits=0 must behave exactly like bits=1 */
        uint32_t coords[3] = { 1, 0, 1 };
        ASSERT_INT_EQ(curve_hilbertnd_xy2d(coords, 3, 0),
                      curve_hilbertnd_xy2d(coords, 3, 1),
                      "hilbertnd bits=0 clamps to 1");
        curve_hilbertnd_d2xy(0, coords, 3, 0);
        x = coords[0]; y = coords[1];
        ASSERT_TRUE(x == 0 && y == 0, "hilbertnd d=0 -> origin");
    }

    {
        uint32_t dummy = 0xDEADu;
        ASSERT_INT_EQ(curve_hilbertnd_xy2d(&dummy, 0, 8), 0,
                      "hilbertnd n_dims=0 -> 0 without reading");
        curve_hilbertnd_d2xy(99, &dummy, 0, 8);
        ASSERT_INT_EQ(dummy, 0xDEADu, "hilbertnd n_dims=0 d2xy does not write");
    }

    /* Out-of-range coordinates are masked onto the grid. */
    ASSERT_INT_EQ(curve_hilbert_xy2d(3, 3, 1), curve_hilbert_xy2d(1, 1, 1),
                  "hilbert2 masks coords to `bits` bits");
    ASSERT_INT_EQ(curve_moore_xy2d(5, 1, 1), curve_moore_xy2d(1, 1, 1),
                  "moore masks coords to `bits` bits");
}

static void test_window_nonfinite(void)
{
    const double w[4] = { 0.0, 0.0, 1.0, 1.0 };
    ASSERT_INT_EQ(curve_z_window(NAN, NAN, w), 0,
                  "NaN window input maps to 0");
    ASSERT_INT_EQ(curve_z_window(-INFINITY, 0.5, w),
                  curve_z_window(0.0, 0.5, w),
                  "-Inf clamps to window min");
    ASSERT_INT_EQ(curve_z_window(INFINITY, 0.5, w),
                  curve_z_window(1.0, 0.5, w),
                  "+Inf clamps to window max");
}

static void test_hilbertnd_bits32(void)
{
    int ok = 1;
    for (int i = 0; i < 20000 && ok; ++i) {
        uint32_t coords[2] = { (uint32_t)i * 2654435761u,
                               (uint32_t)i * 2246822519u + 7u };
        uint64_t d = curve_hilbertnd_xy2d(coords, 2, 32);
        uint32_t back[2];
        curve_hilbertnd_d2xy(d, back, 2, 32);
        if (back[0] != coords[0] || back[1] != coords[1]) ok = 0;
    }
    ASSERT_TRUE(ok, "hilbertnd n=2 bits=32 sampled roundtrip");
}

static void test_batch_matches_scalar(void)
{
    enum { N = 67 }; /* 67 = 16*4 + 3, exercises SIMD body + scalar tail */
    uint32_t x[N], y[N], z[N], d32[N], rx[N], ry[N], rz[N];
    uint64_t d64[N], g[N], gout[N];
    int ok = 1, i;

    for (i = 0; i < N; ++i) {
        x[i] = (uint32_t)i * 2654435761u;
        y[i] = (uint32_t)i * 2246822519u + 13u;
        z[i] = (uint32_t)i * 3266489917u + 29u;
        g[i] = (uint64_t)i * 11400714819323198485ull + 1ull;
    }

    curve_morton2_batch(x, y, d32, 0);
    curve_morton2_batch(x, y, d32, N);
    for (i = 0; i < N && ok; ++i)
        if (d32[i] != curve_morton2(x[i], y[i])) ok = 0;
    ASSERT_TRUE(ok, "morton2_batch matches scalar");

    ok = 1;
    curve_morton2_decode_batch(d32, rx, ry, N);
    for (i = 0; i < N && ok; ++i)
        if (rx[i] != (x[i] & 0xFFFFu) || ry[i] != (y[i] & 0xFFFFu)) ok = 0;
    ASSERT_TRUE(ok, "morton2_decode_batch matches scalar (16-bit axes)");

    ok = 1;
    curve_morton2_64_batch(x, y, d64, N);
    for (i = 0; i < N && ok; ++i)
        if (d64[i] != curve_morton2_64(x[i], y[i])) ok = 0;
    ASSERT_TRUE(ok, "morton2_64_batch matches scalar");

    ok = 1;
    curve_morton2_64_decode_batch(d64, rx, ry, N);
    for (i = 0; i < N && ok; ++i)
        if (rx[i] != x[i] || ry[i] != y[i]) ok = 0;
    ASSERT_TRUE(ok, "morton2_64_decode_batch roundtrip");

    ok = 1;
    curve_morton3_batch(x, y, z, d64, N);
    for (i = 0; i < N && ok; ++i)
        if (d64[i] != curve_morton3(x[i], y[i], z[i])) ok = 0;
    ASSERT_TRUE(ok, "morton3_batch matches scalar");

    ok = 1;
    curve_morton3_decode_batch(d64, rx, ry, rz, N);
    for (i = 0; i < N && ok; ++i)
        if (rx[i] != (x[i] & 0x1FFFFFu) || ry[i] != (y[i] & 0x1FFFFFu) ||
            rz[i] != (z[i] & 0x1FFFFFu))
            ok = 0;
    ASSERT_TRUE(ok, "morton3_decode_batch matches scalar (21-bit axes)");

    ok = 1;
    curve_gray_encode_batch(g, gout, N);
    for (i = 0; i < N && ok; ++i)
        if (gout[i] != curve_gray_encode(g[i])) ok = 0;
    curve_gray_decode_batch(gout, d64, N);
    for (i = 0; i < N && ok; ++i)
        if (d64[i] != g[i]) ok = 0;
    ASSERT_TRUE(ok, "gray batch encode/decode matches scalar");

    ok = 1;
    curve_hilbert_xy2d_batch(x, y, 16, d32, N);
    for (i = 0; i < N && ok; ++i)
        if (d32[i] != curve_hilbert_xy2d(x[i], y[i], 16)) ok = 0;
    curve_hilbert_d2xy_batch(d32, 16, rx, ry, N);
    for (i = 0; i < N && ok; ++i) {
        uint32_t sx, sy;
        curve_hilbert_d2xy(d32[i], 16, &sx, &sy);
        if (rx[i] != sx || ry[i] != sy) ok = 0;
    }
    ASSERT_TRUE(ok, "hilbert2 batch matches scalar");

    ASSERT_TRUE(curve_simd_isa()[0] != '\0', "curve_simd_isa is non-empty");
    ASSERT_INT_EQ(CURVES_VERSION_MAJOR, 1, "version major");
}

int main(void)
{
    RUN_TEST(test_gray_inverse_and_adjacency);
    RUN_TEST(test_morton2_roundtrip);
    RUN_TEST(test_morton3_roundtrip);
    RUN_TEST(test_morton2_64);
    RUN_TEST(test_hilbert2_known_vectors);
    RUN_TEST(test_hilbert2_bijective_all_bits);
    RUN_TEST(test_hilbert2_adjacency);
    RUN_TEST(test_hilbert2_sampled_16bit);
    RUN_TEST(test_hilbert2_64);
    RUN_TEST(test_hilbert3_bijective);
    RUN_TEST(test_hilbertnd_roundtrip);
    RUN_TEST(test_hilbertnd_matches_2d);
    RUN_TEST(test_hilbertnd_prefix_cells);
    RUN_TEST(test_ref_pack_roundtrip);
    RUN_TEST(test_moore_known_path);
    RUN_TEST(test_moore_bijective_continuous);
    RUN_TEST(test_sierpinski_ordering);
    RUN_TEST(test_windowed_matches_discrete);
    RUN_TEST(test_gilbert2_exact_cover_continuous);
    RUN_TEST(test_gilbert3_exact_cover_continuous);
    RUN_TEST(test_clamping);
    RUN_TEST(test_window_nonfinite);
    RUN_TEST(test_hilbertnd_bits32);
    RUN_TEST(test_batch_matches_scalar);
    TEST_SUMMARY();
}
