#include "curves.h"
#include "test_common.h"

#include <limits.h>
#include <stdlib.h>
#include <string.h>

static uint32_t rng_u32(uint32_t *s)
{
    *s = *s * 1664525u + 1013904223u;
    return *s;
}

static int words_equal(const uint64_t *a, const uint64_t *b, size_t n)
{
    size_t i;
    for (i = 0; i < n; i++)
        if (a[i] != b[i]) return 0;
    return 1;
}

static void fill_u64(uint64_t *p, size_t n, uint64_t v)
{
    size_t i;
    for (i = 0; i < n; i++) p[i] = v;
}

static size_t g_sort_words;

static int cmp_keys(const void *pa, const void *pb)
{
    return curve_key_compare((const uint64_t *)pa, (const uint64_t *)pb,
                             g_sort_words);
}

static size_t brute_lower(const uint64_t *keys, size_t n, size_t stride,
                          const uint64_t *query, size_t n_words)
{
    size_t i;
    for (i = 0; i < n; i++) {
        if (curve_key_compare(keys + i * stride, query, n_words) >= 0)
            return i;
    }
    return n;
}

static size_t brute_upper(const uint64_t *keys, size_t n, size_t stride,
                          const uint64_t *query, size_t n_words)
{
    size_t i;
    for (i = 0; i < n; i++) {
        if (curve_key_compare(keys + i * stride, query, n_words) > 0)
            return i;
    }
    return n;
}

static void test_get_bit_and_mask(void)
{
    uint64_t in[2], out[2];
    unsigned bit;
    size_t i;

    in[0] = 0x8000000000000001ull;
    in[1] = 0xC000000000000000ull;
    ASSERT_INT_EQ(curve_key_get_bit(in, 66, 0, &bit), CURVE_OK, "get bit 0");
    ASSERT_INT_EQ(bit, 1u, "MSB set");
    ASSERT_INT_EQ(curve_key_get_bit(in, 66, 63, &bit), CURVE_OK, "get bit 63");
    ASSERT_INT_EQ(bit, 1u, "bit 63 of word0");
    ASSERT_INT_EQ(curve_key_get_bit(in, 66, 64, &bit), CURVE_OK, "get bit 64");
    ASSERT_INT_EQ(bit, 1u, "MSB of word1");
    ASSERT_INT_EQ(curve_key_get_bit(in, 66, 65, &bit), CURVE_OK, "get bit 65");
    ASSERT_INT_EQ(bit, 1u, "second bit of word1");
    ASSERT_INT_EQ(curve_key_get_bit(in, 66, 66, &bit), CURVE_EINVAL,
                  "past end");
    ASSERT_INT_EQ(curve_key_get_bit(NULL, 66, 0, &bit), CURVE_EINVAL, "null");

    fill_u64(out, 2, 0xDEADBEEFDEADBEEFull);
    ASSERT_INT_EQ(curve_key_prefix_mask(in, 66, 0, out, 2), CURVE_OK,
                  "mask 0");
    ASSERT_TRUE(out[0] == 0 && out[1] == 0, "empty prefix zeros");

    ASSERT_INT_EQ(curve_key_prefix_mask(in, 66, 64, out, 2), CURVE_OK,
                  "mask 64");
    ASSERT_INT_EQ(out[0], in[0], "first word kept");
    ASSERT_INT_EQ(out[1], 0ull, "rest cleared, no 64-bit shift");

    ASSERT_INT_EQ(curve_key_prefix_mask(in, 66, 65, out, 2), CURVE_OK,
                  "mask 65");
    ASSERT_INT_EQ(out[0], in[0], "word0 kept");
    ASSERT_INT_EQ(out[1], 0x8000000000000000ull, "one bit of word1");

    ASSERT_INT_EQ(curve_key_prefix_mask(in, 66, 66, out, 2), CURVE_OK,
                  "mask all");
    ASSERT_TRUE(words_equal(in, out, 2), "full prefix is a copy");

    ASSERT_INT_EQ(curve_key_prefix_mask(in, 66, 67, out, 2), CURVE_EINVAL,
                  "prefix too long");
    ASSERT_INT_EQ(curve_key_prefix_mask(in, 66, 8, out, 1), CURVE_ENOSPACE,
                  "short out");

    /* Unused tail of a non-multiple-of-64 key is cleared. */
    in[0] = 0xFF000000000000FFull;
    ASSERT_INT_EQ(curve_key_prefix_mask(in, 8, 8, out, 1), CURVE_OK,
                  "K=8 full");
    ASSERT_INT_EQ(out[0], 0xFF00000000000000ull, "unused low bits zeroed");

    for (i = 0; i < 64; i++) {
        in[0] = 1ull << (63u - (unsigned)i);
        curve_key_get_bit(in, 64, i, &bit);
        ASSERT_INT_EQ(bit, 1u, "each bit position");
    }
}

static void test_common_prefix(void)
{
    uint64_t a[2], b[2];
    size_t p;

    a[0] = b[0] = 0x0123456789ABCDEFull;
    a[1] = b[1] = 0xF000000000000000ull;
    ASSERT_INT_EQ(curve_key_common_prefix_bits(a, b, 68, &p), CURVE_OK,
                  "identical");
    ASSERT_INT_EQ(p, 68, "full 68");

    b[1] = 0x7000000000000000ull; /* bit 68's word: high nibble 0111 vs 1111 */
    ASSERT_INT_EQ(curve_key_common_prefix_bits(a, b, 68, &p), CURVE_OK,
                  "differ at bit 64");
    ASSERT_INT_EQ(p, 64, "first 64 match");

    b[0] = a[0] ^ (1ull << 3);
    b[1] = a[1];
    ASSERT_INT_EQ(curve_key_common_prefix_bits(a, b, 68, &p), CURVE_OK,
                  "differ in word0");
    ASSERT_INT_EQ(p, 60, "clz of bit 3 from LSB of word = bit 60 from MSB");

    a[0] = 0; b[0] = ~0ull;
    ASSERT_INT_EQ(curve_key_common_prefix_bits(a, b, 64, &p), CURVE_OK,
                  "total mismatch");
    ASSERT_INT_EQ(p, 0, "no shared bits");

    /* Unused tail bits must not affect the prefix length. */
    a[0] = b[0] = 0xAB00000000000000ull;
    a[0] |= 0x1; /* junk in unused low bits of a K=8 key */
    ASSERT_INT_EQ(curve_key_common_prefix_bits(a, b, 8, &p), CURVE_OK,
                  "ignore unused");
    ASSERT_INT_EQ(p, 8, "unused xor ignored");

    ASSERT_INT_EQ(curve_key_common_prefix_bits(NULL, b, 8, &p), CURVE_EINVAL,
                  "null");
}

static void test_matches_successor_interval(void)
{
    uint64_t key[2], pref[2], lower[2], upper[2], sent;
    int match, unbounded;

    key[0] = 0xA000000000000000ull;
    key[1] = 0;
    curve_key_prefix_mask(key, 64, 4, pref, 2);
    ASSERT_INT_EQ(curve_key_prefix_matches(key, pref, 64, 4, &match),
                  CURVE_OK, "matches");
    ASSERT_INT_EQ(match, 1, "4-bit prefix hits");
    pref[0] = 0xB000000000000000ull;
    ASSERT_INT_EQ(curve_key_prefix_matches(key, pref, 64, 4, &match),
                  CURVE_OK, "no match");
    ASSERT_INT_EQ(match, 0, "different nibble");
    ASSERT_INT_EQ(curve_key_prefix_matches(key, pref, 64, 0, &match),
                  CURVE_OK, "empty prefix");
    ASSERT_INT_EQ(match, 1, "empty matches everything");

    sent = 0xDEADBEEFDEADBEEFull;
    key[1] = 0x1234000000000000ull;
    fill_u64(upper, 2, sent);
    ASSERT_INT_EQ(curve_key_prefix_successor(key, 128, 4, upper, 2), CURVE_OK,
                  "succ of 0xA");
    ASSERT_INT_EQ(upper[0], 0xB000000000000000ull, "0xA + 1 = 0xB");
    ASSERT_INT_EQ(upper[1], 0ull, "rest of 128-bit key zero");
    key[1] = 0;

    /* All-ones prefix: ERANGE and `out` left unchanged. */
    key[0] = ~0ull;
    fill_u64(upper, 2, sent);
    ASSERT_INT_EQ(curve_key_prefix_successor(key, 64, 64, upper, 2),
                  CURVE_ERANGE, "all-ones 64");
    ASSERT_INT_EQ(upper[0], sent, "ERANGE does not write");

    key[0] = ~0ull;
    key[1] = 0xFC00000000000000ull; /* 6 high bits set */
    fill_u64(upper, 2, sent);
    ASSERT_INT_EQ(curve_key_prefix_successor(key, 70, 70, upper, 2),
                  CURVE_ERANGE, "all-ones 70");
    ASSERT_TRUE(upper[0] == sent && upper[1] == sent, "70-bit ERANGE no write");

    ASSERT_INT_EQ(curve_key_prefix_successor(key, 70, 0, upper, 2),
                  CURVE_ERANGE, "empty prefix no successor");

    /* Carry across a word boundary: prefix 64 of 0xFF..FE. */
    key[0] = ~0ull - 1ull;
    key[1] = 0;
    ASSERT_INT_EQ(curve_key_prefix_successor(key, 70, 64, upper, 2), CURVE_OK,
                  "carry to next word");
    ASSERT_TRUE(upper[0] == ~0ull && upper[1] == 0ull, "64-bit wrap to 0xFF..FF");

    /* Non-terminal prefix 0b10. */
    key[0] = 0x8000000000000000ull;
    key[1] = 0;
    ASSERT_INT_EQ(curve_key_prefix_interval(key, 64, 2, lower, upper, 2,
                                            &unbounded),
                  CURVE_OK, "interval 0b10");
    ASSERT_INT_EQ(unbounded, 0, "bounded");
    ASSERT_INT_EQ(lower[0], 0x8000000000000000ull, "lower 0b10");
    ASSERT_INT_EQ(upper[0], 0xC000000000000000ull, "upper 0b11");

    /* Terminal prefix 0b11: successor carries out. */
    key[0] = 0xC000000000000000ull;
    ASSERT_INT_EQ(curve_key_prefix_interval(key, 64, 2, lower, upper, 2,
                                            &unbounded),
                  CURVE_OK, "terminal 2-bit prefix");
    ASSERT_INT_EQ(unbounded, 1, "0b11 is terminal");
    ASSERT_INT_EQ(lower[0], 0xC000000000000000ull, "lower kept");
    ASSERT_INT_EQ(upper[0], 0ull, "upper zeroed when unbounded");

    key[0] = 0;
    ASSERT_INT_EQ(curve_key_prefix_interval(key, 64, 0, lower, upper, 2,
                                            &unbounded),
                  CURVE_OK, "whole domain");
    ASSERT_INT_EQ(unbounded, 1, "empty prefix unbounded");
    ASSERT_INT_EQ(lower[0], 0ull, "lower zero");

    ASSERT_INT_EQ(curve_key_prefix_interval(key, 64, 1, lower, NULL, 2,
                                            &unbounded),
                  CURVE_EINVAL, "null upper");

    ASSERT_INT_EQ(curve_key_prefix_is_cell(4, 2, 8), 1, "complete level");
    ASSERT_INT_EQ(curve_key_prefix_is_cell(3, 2, 8), 0, "partial");
    ASSERT_INT_EQ(curve_key_prefix_is_cell(0, 2, 8), 1, "root is a cell");
    ASSERT_INT_EQ(curve_key_cell_level(10, 2), 5, "level");
    ASSERT_INT_EQ(curve_key_cell_prefix_bits(5, 2), 10, "prefix bits");
}

static void test_binary_search_vs_brute(void)
{
    enum { N = 200, Q = 400, STRIDE = 3 };
    uint64_t keys[N * STRIDE];
    uint64_t packed[N];
    uint64_t query[2];
    uint32_t seed = 1;
    size_t i, idx, brute, n_words;
    int ok;

    /* One-word keys, stride 1. */
    n_words = 1;
    for (i = 0; i < N; i++)
        packed[i] = ((uint64_t)rng_u32(&seed) << 32) | rng_u32(&seed);
    g_sort_words = n_words;
    qsort(packed, N, sizeof(uint64_t), cmp_keys);

    ok = 1;
    for (i = 0; i < Q && ok; i++) {
        query[0] = ((uint64_t)rng_u32(&seed) << 32) | rng_u32(&seed);
        if (i == 0) query[0] = packed[0];
        if (i == 1) query[0] = packed[N - 1];
        if (i == 2) query[0] = 0;
        if (i == 3) query[0] = ~0ull;
        curve_key_lower_bound(packed, N, 1, query, 1, &idx);
        brute = brute_lower(packed, N, 1, query, 1);
        if (idx != brute) ok = 0;
        curve_key_upper_bound(packed, N, 1, query, 1, &idx);
        brute = brute_upper(packed, N, 1, query, 1);
        if (idx != brute) ok = 0;
    }
    ASSERT_TRUE(ok, "1-word lower/upper vs brute");

    /* Two-word keys packed with stride 3. */
    n_words = 2;
    for (i = 0; i < N; i++) {
        keys[i * STRIDE + 0] = ((uint64_t)rng_u32(&seed) << 32) | rng_u32(&seed);
        keys[i * STRIDE + 1] = ((uint64_t)rng_u32(&seed) << 32) | rng_u32(&seed);
        keys[i * STRIDE + 2] = 0xDEADBEEFDEADBEEFull; /* padding */
    }
    g_sort_words = n_words;
    qsort(keys, N, STRIDE * sizeof(uint64_t), cmp_keys);

    ok = 1;
    for (i = 0; i < Q && ok; i++) {
        query[0] = ((uint64_t)rng_u32(&seed) << 32) | rng_u32(&seed);
        query[1] = ((uint64_t)rng_u32(&seed) << 32) | rng_u32(&seed);
        curve_key_lower_bound(keys, N, STRIDE, query, 2, &idx);
        brute = brute_lower(keys, N, STRIDE, query, 2);
        if (idx != brute) ok = 0;
        curve_key_upper_bound(keys, N, STRIDE, query, 2, &idx);
        brute = brute_upper(keys, N, STRIDE, query, 2);
        if (idx != brute) ok = 0;
    }
    ASSERT_TRUE(ok, "2-word stride-3 lower/upper vs brute");

    {
        size_t z;
        ASSERT_INT_EQ(curve_key_lower_bound(NULL, 0, 1, NULL, 1, &z),
                      CURVE_OK, "empty array");
        ASSERT_INT_EQ(z, 0, "empty -> 0");
        ASSERT_INT_EQ(curve_key_lower_bound(packed, N, 1, query, 0, &z),
                      CURVE_EINVAL, "zero n_words");
        ASSERT_INT_EQ(curve_key_lower_bound(packed, N, 0, query, 1, &z),
                      CURVE_EINVAL, "stride < n_words");
    }
}

static void test_prefix_range_partition(void)
{
    enum { N = 256 };
    uint64_t keys[N], scratch[2];
    uint32_t coords[2], x, y;
    size_t i, begin, end, pbits;
    int ok, match;

    for (y = 0; y < 16; y++) {
        for (x = 0; x < 16; x++) {
            coords[0] = x;
            coords[1] = y;
            curve_mortonnd_encode_words(coords, 2, 4, &keys[y * 16 + x], 1);
        }
    }
    g_sort_words = 1;
    qsort(keys, N, sizeof(uint64_t), cmp_keys);

    ok = 1;
    for (pbits = 0; pbits <= 8 && ok; pbits++) {
        for (i = 0; i < N && ok; i++) {
            size_t j, brute_lo, brute_hi;
            curve_key_prefix_range(keys, N, 1, &keys[i], 8, pbits,
                                   scratch, 2, &begin, &end);
            brute_lo = N;
            brute_hi = 0;
            for (j = 0; j < N; j++) {
                curve_key_prefix_matches(&keys[j], &keys[i], 8, pbits, &match);
                if (match) {
                    if (j < brute_lo) brute_lo = j;
                    brute_hi = j + 1;
                }
            }
            if (pbits == 0) {
                if (begin != 0 || end != N) ok = 0;
            } else if (brute_lo == N) {
                if (begin != end) ok = 0;
            } else if (begin != brute_lo || end != brute_hi) {
                ok = 0;
            }
            /* Matching keys are a contiguous slice of the sorted array. */
            for (j = begin; j < end && ok; j++) {
                curve_key_prefix_matches(&keys[j], &keys[i], 8, pbits, &match);
                if (!match) ok = 0;
            }
        }
    }
    ASSERT_TRUE(ok, "prefix ranges partition the sorted domain");

    /* All-ones key: prefix_range runs to n_keys (unbounded). */
    {
        uint64_t ones = ~0ull;
        uint64_t arr[3];
        arr[0] = 0;
        arr[1] = 0x8000000000000000ull;
        arr[2] = ~0ull;
        curve_key_prefix_range(arr, 3, 1, &ones, 64, 64, scratch, 2,
                               &begin, &end);
        ASSERT_INT_EQ(begin, 2, "all-ones lower");
        ASSERT_INT_EQ(end, 3, "unbounded end = n_keys");
    }

    ASSERT_INT_EQ(curve_key_prefix_range(keys, N, 1, keys, 8, 4, scratch, 1,
                                         &begin, &end),
                  CURVE_ENOSPACE, "scratch too small");
}

static void test_cell_init_nav(void)
{
    curve_cell cell, parent, child, sib;
    uint64_t key[1], pwords[1], cwords[1], swords[1];
    uint32_t coords[2] = { 3u, 5u };
    uint64_t idx;
    size_t nchild;
    size_t i;
    curve_status st;

    ASSERT_INT_EQ(curve_cell_init(&cell, 2, 4, 3), CURVE_EINVAL,
                  "partial level");
    ASSERT_INT_EQ(curve_cell_init(&cell, 2, 4, 10), CURVE_EINVAL,
                  "prefix past K");
    ASSERT_INT_EQ(curve_cell_init(&cell, 2, 4, 0), CURVE_OK, "root");
    ASSERT_INT_EQ(curve_cell_depth(&cell), 0, "root depth");
    ASSERT_INT_EQ(curve_cell_is_leaf(&cell), 0, "root not leaf");
    ASSERT_INT_EQ(curve_cell_child_count(&cell, &nchild), CURVE_OK,
                  "root fanout");
    ASSERT_INT_EQ(nchild, 4, "2D fanout 4");

    curve_mortonnd_encode_words(coords, 2, 4, key, 1);
    ASSERT_INT_EQ(curve_cell_init(&cell, 2, 4, 4), CURVE_OK, "depth 2");
    st = curve_cell_parent(&cell, key, &parent, pwords, 1);
    ASSERT_INT_EQ(st, CURVE_OK, "parent");
    ASSERT_INT_EQ(parent.prefix_bits, 2, "parent is depth 1");

    st = curve_cell_child_index(&cell, key, &idx);
    ASSERT_INT_EQ(st, CURVE_OK, "child index");
    st = curve_cell_child(&parent, pwords, idx, &child, cwords, 1);
    ASSERT_INT_EQ(st, CURVE_OK, "re-child");
    ASSERT_INT_EQ(child.prefix_bits, 4, "child depth restored");
    {
        uint64_t masked[1], masked2[1];
        curve_key_prefix_mask(key, 8, 4, masked, 1);
        curve_key_prefix_mask(cwords, 8, 4, masked2, 1);
        ASSERT_TRUE(masked[0] == masked2[0], "child prefix matches");
    }

    /* Walk all 4 children of the root; they cover child indices 0..3. */
    ASSERT_INT_EQ(curve_cell_init(&cell, 2, 4, 0), CURVE_OK, "root again");
    key[0] = 0;
    for (i = 0; i < 4; i++) {
        st = curve_cell_child(&cell, key, i, &child, cwords, 1);
        ASSERT_INT_EQ(st, CURVE_OK, "child i");
        st = curve_cell_child_index(&child, cwords, &idx);
        ASSERT_INT_EQ(idx, (uint64_t)i, "index roundtrip");
    }
    ASSERT_INT_EQ(curve_cell_child(&cell, key, 4, &child, cwords, 1),
                  CURVE_EINVAL, "child 4 invalid");

    /* Siblings. */
    ASSERT_INT_EQ(curve_cell_child(&cell, key, 1, &child, cwords, 1),
                  CURVE_OK, "child 1");
    ASSERT_INT_EQ(curve_cell_sibling(&child, cwords, 1, &sib, swords, 1),
                  CURVE_OK, "next sib");
    curve_cell_child_index(&sib, swords, &idx);
    ASSERT_INT_EQ(idx, 2, "sibling +1 -> index 2");
    ASSERT_INT_EQ(curve_cell_sibling(&child, cwords, -1, &sib, swords, 1),
                  CURVE_OK, "prev sib");
    curve_cell_child_index(&sib, swords, &idx);
    ASSERT_INT_EQ(idx, 0, "sibling -1 -> index 0");
    ASSERT_INT_EQ(curve_cell_sibling(&child, cwords, -2, &sib, swords, 1),
                  CURVE_ERANGE, "past first sibling");
    ASSERT_INT_EQ(curve_cell_sibling(&child, cwords, 3, &sib, swords, 1),
                  CURVE_ERANGE, "past last sibling");
    ASSERT_INT_EQ(curve_cell_sibling(&child, cwords, 0, &sib, swords, 1),
                  CURVE_OK, "delta 0 identity");
    ASSERT_INT_EQ(sib.prefix_bits, child.prefix_bits, "identity cell");
    ASSERT_INT_EQ(curve_cell_sibling(&child, cwords, INT_MIN, &sib, swords, 1),
                  CURVE_ERANGE, "INT_MIN");

    /* In-place sibling (sib_prefix aliases prefix). */
    ASSERT_INT_EQ(curve_cell_sibling(&child, cwords, 1, &sib, cwords, 1),
                  CURVE_OK, "in-place sibling");
    curve_cell_child_index(&sib, cwords, &idx);
    ASSERT_INT_EQ(idx, 2, "in-place +1");

    /* Root has no parent; leaf has no child. */
    ASSERT_INT_EQ(curve_cell_init(&cell, 2, 4, 0), CURVE_OK, "root");
    key[0] = 0;
    ASSERT_INT_EQ(curve_cell_parent(&cell, key, &parent, pwords, 1),
                  CURVE_ERANGE, "root parent");
    ASSERT_INT_EQ(curve_cell_init(&cell, 2, 4, 8), CURVE_OK, "leaf");
    ASSERT_TRUE(curve_cell_is_leaf(&cell), "is leaf");
    ASSERT_INT_EQ(curve_cell_child_count(&cell, &nchild), CURVE_OK, "leaf count");
    ASSERT_INT_EQ(nchild, 0, "leaf has 0 children");
    ASSERT_INT_EQ(curve_cell_child(&cell, key, 0, &child, cwords, 1),
                  CURVE_ERANGE, "leaf child");

    {
        uint64_t wide[5];
        size_t big;
        fill_u64(wide, 5, 0);
        ASSERT_INT_EQ(curve_cell_init(&cell, 65, 4, 0), CURVE_OK, "65D root");
        ASSERT_INT_EQ(curve_cell_child_count(&cell, &big), CURVE_EUNSUPPORTED,
                      "65D child_count");
        ASSERT_INT_EQ(curve_cell_child(&cell, wide, 0, &child, wide, 5),
                      CURVE_EUNSUPPORTED, "65D child");
    }
}

static int aabb_contains(const uint32_t *mins, const uint32_t *maxs,
                         const uint32_t *c, size_t n)
{
    size_t i;
    for (i = 0; i < n; i++)
        if (c[i] < mins[i] || c[i] > maxs[i]) return 0;
    return 1;
}

static int test_grid_cells(int hilbert, size_t n, unsigned bits)
{
    uint32_t coords[8], other[8], scratch[8];
    uint32_t mins[8], maxs[8], pmins[8], pmaxs[8];
    uint64_t key[1], pref[1], pwords[1], cwords[1];
    curve_cell cell, parent, child;
    size_t side, total, j, d, level, max_level;
    uint32_t t;
    int match;
    curve_status st;

    if (n > 8u) return 0;
    side = (size_t)1u << bits;
    total = 1;
    for (d = 0; d < n; d++) total *= side;
    if (total > 4096u) return 0;
    max_level = (size_t)bits;

    for (t = 0; t < (uint32_t)total; t++) {
        uint32_t v = t;
        for (d = 0; d < n; d++) {
            coords[d] = v % (uint32_t)side;
            v /= (uint32_t)side;
        }
        if (hilbert)
            curve_hilbertnd_encode_words(coords, n, bits, key, 1, scratch, 8);
        else
            curve_mortonnd_encode_words(coords, n, bits, key, 1);

        for (level = 0; level <= max_level; level++) {
            size_t pbits = level * n;
            uint32_t side_cell = (level >= (size_t)bits)
                ? 1u : (1u << (bits - (unsigned)level));
            st = curve_cell_init(&cell, n, bits, pbits);
            if (st != CURVE_OK) return 0;
            curve_key_prefix_mask(key, n * (size_t)bits, pbits, pref, 1);
            if (hilbert)
                st = curve_hilbert_cell_bounds(&cell, pref, mins, maxs, 8,
                                               scratch, 8);
            else
                st = curve_morton_cell_bounds(&cell, pref, mins, maxs, 8);
            if (st != CURVE_OK) return 0;
            if (!aabb_contains(mins, maxs, coords, n)) return 0;
            for (d = 0; d < n; d++)
                if (maxs[d] - mins[d] + 1u != side_cell) return 0;

            /* Every grid point in the AABB shares this prefix; none outside. */
            for (j = 0; j < total; j++) {
                uint32_t w = (uint32_t)j;
                uint64_t okw[1];
                for (d = 0; d < n; d++) {
                    other[d] = w % (uint32_t)side;
                    w /= (uint32_t)side;
                }
                if (hilbert)
                    curve_hilbertnd_encode_words(other, n, bits, okw, 1,
                                                 scratch, 8);
                else
                    curve_mortonnd_encode_words(other, n, bits, okw, 1);
                curve_key_prefix_matches(okw, pref, n * (size_t)bits, pbits,
                                         &match);
                if (aabb_contains(mins, maxs, other, n)) {
                    if (!match) return 0;
                } else if (match) {
                    return 0;
                }
            }

            if (level == 0) continue;
            st = curve_cell_parent(&cell, pref, &parent, pwords, 1);
            if (st != CURVE_OK) return 0;
            if (hilbert)
                curve_hilbert_cell_bounds(&parent, pwords, pmins, pmaxs, 8,
                                          scratch, 8);
            else
                curve_morton_cell_bounds(&parent, pwords, pmins, pmaxs, 8);
            if (!aabb_contains(pmins, pmaxs, mins, n)) return 0;
            if (!aabb_contains(pmins, pmaxs, maxs, n)) return 0;

            if (n <= 3u && !curve_cell_is_leaf(&parent)) {
                uint64_t cidx;
                size_t nch = 0, vol = 1, cvol = 0;
                for (d = 0; d < n; d++)
                    vol *= (size_t)(pmaxs[d] - pmins[d] + 1u);
                curve_cell_child_count(&parent, &nch);
                for (cidx = 0; cidx < (uint64_t)nch; cidx++) {
                    uint32_t cmins[8], cmaxs[8];
                    curve_cell_child(&parent, pwords, cidx, &child, cwords, 1);
                    if (hilbert)
                        curve_hilbert_cell_bounds(&child, cwords, cmins, cmaxs,
                                                   8, scratch, 8);
                    else
                        curve_morton_cell_bounds(&child, cwords, cmins, cmaxs, 8);
                    if (!aabb_contains(pmins, pmaxs, cmins, n)) return 0;
                    if (!aabb_contains(pmins, pmaxs, cmaxs, n)) return 0;
                    {
                        size_t cv = 1;
                        for (d = 0; d < n; d++)
                            cv *= (size_t)(cmaxs[d] - cmins[d] + 1u);
                        cvol += cv;
                    }
                }
                if (cvol != vol) return 0;
            }
        }
    }
    return 1;
}

static void test_morton_cell_bounds_grid(void)
{
    ASSERT_TRUE(test_grid_cells(0, 2, 1), "morton 2D bits=1");
    ASSERT_TRUE(test_grid_cells(0, 2, 3), "morton 2D bits=3");
    ASSERT_TRUE(test_grid_cells(0, 2, 4), "morton 2D bits=4");
    ASSERT_TRUE(test_grid_cells(0, 3, 2), "morton 3D bits=2");
}

static void test_hilbert_cell_bounds_grid(void)
{
    ASSERT_TRUE(test_grid_cells(1, 2, 1), "hilbert 2D bits=1");
    ASSERT_TRUE(test_grid_cells(1, 2, 3), "hilbert 2D bits=3");
    ASSERT_TRUE(test_grid_cells(1, 2, 4), "hilbert 2D bits=4");
    ASSERT_TRUE(test_grid_cells(1, 3, 2), "hilbert 3D bits=2");
}

static void test_wide_prefix_and_root_bounds(void)
{
    uint32_t coords[256], scratch[256], mins[256], maxs[256];
    uint64_t morton[16], hilbert[16], masked[16];
    curve_cell cell;
    size_t i, p;
    int match;

    for (i = 0; i < 256; i++) coords[i] = (uint32_t)(i & 15u);
    ASSERT_INT_EQ(curve_mortonnd_encode_words(coords, 256, 4, morton, 16),
                  CURVE_OK, "256x4 morton");
    ASSERT_INT_EQ(curve_hilbertnd_encode_words(coords, 256, 4, hilbert, 16,
                                               scratch, 256),
                  CURVE_OK, "256x4 hilbert");
    ASSERT_INT_EQ(curve_key_common_prefix_bits(morton, morton, 1024, &p),
                  CURVE_OK, "self prefix");
    ASSERT_INT_EQ(p, 1024, "1024 bits");

    ASSERT_INT_EQ(curve_cell_init(&cell, 256, 4, 0), CURVE_OK, "256D root");
    ASSERT_INT_EQ(curve_morton_cell_bounds(&cell, morton, mins, maxs, 256),
                  CURVE_OK, "morton root bounds");
    ASSERT_TRUE(mins[0] == 0 && maxs[0] == 15 && mins[255] == 0 &&
                maxs[255] == 15,
                "root is the full 4-bit domain");
    ASSERT_INT_EQ(curve_hilbert_cell_bounds(&cell, hilbert, mins, maxs, 256,
                                            scratch, 256),
                  CURVE_OK, "hilbert root bounds");
    ASSERT_TRUE(mins[0] == 0 && maxs[0] == 15, "hilbert root full domain");

    /* Depth-1 cell: 256 prefix bits = 4 words. */
    ASSERT_INT_EQ(curve_cell_init(&cell, 256, 4, 256), CURVE_OK, "depth 1");
    curve_key_prefix_mask(morton, 1024, 256, masked, 16);
    ASSERT_INT_EQ(curve_morton_cell_bounds(&cell, masked, mins, maxs, 256),
                  CURVE_OK, "morton depth-1");
    ASSERT_TRUE(aabb_contains(mins, maxs, coords, 256),
                "point inside morton depth-1 cell");
    curve_key_prefix_matches(morton, masked, 1024, 256, &match);
    ASSERT_INT_EQ(match, 1, "key matches its depth-1 prefix");

    curve_key_prefix_mask(hilbert, 1024, 256, masked, 16);
    ASSERT_INT_EQ(curve_hilbert_cell_bounds(&cell, masked, mins, maxs, 256,
                                            scratch, 256),
                  CURVE_OK, "hilbert depth-1, W=16 > 32-word old cap");
    ASSERT_TRUE(aabb_contains(mins, maxs, coords, 256),
                "point inside hilbert depth-1 cell");
}

static void test_prefix_failures(void)
{
    uint64_t w[2] = { 0, 0 };
    uint64_t out[2];
    size_t p;
    unsigned bit;
    int match;
    curve_cell cell;

    ASSERT_INT_EQ(curve_key_get_bit(w, 0, 0, &bit), CURVE_EINVAL, "K=0");
    ASSERT_INT_EQ(curve_key_common_prefix_bits(w, w, 0, &p), CURVE_EINVAL,
                  "common K=0");
    ASSERT_INT_EQ(curve_key_prefix_mask(w, 8, 9, out, 2), CURVE_EINVAL,
                  "mask too long");
    ASSERT_INT_EQ(curve_key_prefix_matches(w, w, 8, 9, &match), CURVE_EINVAL,
                  "matches too long");
    ASSERT_INT_EQ(curve_key_prefix_interval(w, 8, 1, out, out, 2, NULL),
                  CURVE_EINVAL, "null unbounded");
    ASSERT_INT_EQ(curve_cell_init(NULL, 2, 4, 0), CURVE_EINVAL, "null cell");
    ASSERT_INT_EQ(curve_cell_init(&cell, 0, 4, 0), CURVE_EINVAL, "zero dims");
}

int main(void)
{
    RUN_TEST(test_get_bit_and_mask);
    RUN_TEST(test_common_prefix);
    RUN_TEST(test_matches_successor_interval);
    RUN_TEST(test_binary_search_vs_brute);
    RUN_TEST(test_prefix_range_partition);
    RUN_TEST(test_cell_init_nav);
    RUN_TEST(test_morton_cell_bounds_grid);
    RUN_TEST(test_hilbert_cell_bounds_grid);
    RUN_TEST(test_wide_prefix_and_root_bounds);
    RUN_TEST(test_prefix_failures);
    TEST_SUMMARY();
}
