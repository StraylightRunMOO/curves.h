/* Phase G: Gray-coded Morton, Peano 2D, mixed-bit Morton. */
#include "curves.h"
#include "test_common.h"

#include <stdlib.h>
#include <string.h>

static uint32_t gray32(uint32_t x)
{
    return x ^ (x >> 1);
}

static void test_gmorton2_matches_gray_then_morton(void)
{
    uint32_t x, y, d, rx, ry;
    for (x = 0; x < 64; x++) {
        for (y = 0; y < 64; y++) {
            d = curve_gmorton2(x, y);
            ASSERT_INT_EQ(d, curve_morton2(gray32(x), gray32(y)),
                          "gmorton2 = morton2(gray x, gray y)");
            curve_gmorton2_decode(d, &rx, &ry);
            ASSERT_INT_EQ(rx, x, "gmorton2 decode x");
            ASSERT_INT_EQ(ry, y, "gmorton2 decode y");
        }
    }
    ASSERT_INT_EQ(curve_gmorton2(0, 0), 0, "origin");
    ASSERT_INT_EQ(curve_gmorton2(0x1FFFFu, 0), curve_gmorton2(0xFFFFu, 0),
                  "gmorton2 masks to 16 bits");
}

static void test_gmorton2_64_roundtrip(void)
{
    uint32_t i, x, y, rx, ry;
    uint64_t d;
    for (i = 0; i < 4096; i++) {
        x = i * 2654435761u;
        y = i * 2246822519u + 13u;
        d = curve_gmorton2_64(x, y);
        ASSERT_INT_EQ(d, curve_morton2_64(gray32(x), gray32(y)),
                      "gmorton2_64 = morton2_64(gray, gray)");
        curve_gmorton2_64_decode(d, &rx, &ry);
        ASSERT_INT_EQ(rx, x, "gmorton2_64 decode x");
        ASSERT_INT_EQ(ry, y, "gmorton2_64 decode y");
    }
}

static void test_gmorton3_roundtrip(void)
{
    uint32_t x, y, z, rx, ry, rz;
    uint64_t d;
    for (x = 0; x < 16; x++) {
        for (y = 0; y < 16; y++) {
            for (z = 0; z < 16; z++) {
                d = curve_gmorton3(x, y, z);
                ASSERT_INT_EQ(d, curve_morton3(gray32(x), gray32(y), gray32(z)),
                              "gmorton3 = morton3(gray)");
                curve_gmorton3_decode(d, &rx, &ry, &rz);
                ASSERT_INT_EQ(rx, x, "gmorton3 x");
                ASSERT_INT_EQ(ry, y, "gmorton3 y");
                ASSERT_INT_EQ(rz, z, "gmorton3 z");
            }
        }
    }
    ASSERT_INT_EQ(curve_gmorton3(0x3FFFFFu, 0, 0),
                  curve_gmorton3(0x1FFFFFu, 0, 0),
                  "gmorton3 masks to 21 bits");
}

static void test_gmortonnd_wide(void)
{
    uint32_t coords[4] = { 7u, 13u, 1u, 9u }, back[4];
    uint64_t words[1];
    uint32_t gcoords[4];
    size_t i;
    ASSERT_INT_EQ(curve_gmortonnd_encode_words(coords, 4, 8, words, 1),
                  CURVE_OK, "gmortonnd encode");
    for (i = 0; i < 4; i++) gcoords[i] = gray32(coords[i] & 0xFFu);
    {
        uint64_t expect[1];
        ASSERT_INT_EQ(curve_mortonnd_encode_words(gcoords, 4, 8, expect, 1),
                      CURVE_OK, "morton of gray coords");
        ASSERT_INT_EQ(words[0], expect[0], "gmortonnd word0");
    }
    ASSERT_INT_EQ(curve_gmortonnd_decode_words(words, 1, 4, 8, back, 4),
                  CURVE_OK, "gmortonnd decode");
    for (i = 0; i < 4; i++)
        ASSERT_INT_EQ(back[i], coords[i], "gmortonnd roundtrip axis");
    {
        uint32_t c2[3] = { 0x00ABCDEFu, 0x00123456u, 0x000000FFu }, b2[3];
        uint64_t w2[2], e2[2];
        uint32_t g2[3];
        ASSERT_INT_EQ(curve_gmortonnd_encode_words(c2, 3, 32, w2, 2),
                      CURVE_OK, "gmortonnd 3x32 encode");
        g2[0] = gray32(c2[0]);
        g2[1] = gray32(c2[1]);
        g2[2] = gray32(c2[2]);
        ASSERT_INT_EQ(curve_mortonnd_encode_words(g2, 3, 32, e2, 2),
                      CURVE_OK, "morton 3x32 of gray");
        ASSERT_INT_EQ(w2[0], e2[0], "gmortonnd 3x32 word0");
        ASSERT_INT_EQ(w2[1], e2[1], "gmortonnd 3x32 word1");
        ASSERT_INT_EQ(curve_gmortonnd_decode_words(w2, 2, 3, 32, b2, 3),
                      CURVE_OK, "gmortonnd 3x32 decode");
        ASSERT_INT_EQ(b2[0], c2[0], "gmortonnd 3x32 axis0");
        ASSERT_INT_EQ(b2[1], c2[1], "gmortonnd 3x32 axis1");
        ASSERT_INT_EQ(b2[2], c2[2], "gmortonnd 3x32 axis2");
    }
    ASSERT_INT_EQ(curve_gmortonnd_encode_words(NULL, 2, 8, words, 1),
                  CURVE_EINVAL, "null coords");
    ASSERT_INT_EQ(curve_gmortonnd_encode_words(coords, 4, 8, words, 0),
                  CURVE_ENOSPACE, "short key");
}

static void test_peano_order1(void)
{
    /* Column snake, origin BL, x right, y up. */
    static const uint32_t expect[3][3] = {
        { 0, 1, 2 },
        { 5, 4, 3 },
        { 6, 7, 8 }
    };
    uint32_t x, y, rx, ry, d;
    for (x = 0; x < 3; x++) {
        for (y = 0; y < 3; y++) {
            d = curve_peano_xy2d(x, y, 1);
            ASSERT_INT_EQ(d, expect[x][y], "peano order-1 generator");
            curve_peano_d2xy(d, 1, &rx, &ry);
            ASSERT_INT_EQ(rx, x, "peano1 decode x");
            ASSERT_INT_EQ(ry, y, "peano1 decode y");
        }
    }
}

static void test_peano_bijection_adjacency(void)
{
    int digits;
    for (digits = 1; digits <= 5; digits++) {
        uint32_t n = 1, side, N, i, x, y, d, rx, ry, px = 0, py = 0;
        uint8_t *seen;
        int k;
        for (k = 0; k < digits; k++) n *= 3u;
        side = n;
        N = side * side;
        seen = (uint8_t *)calloc(N, 1);
        ASSERT_TRUE(seen != NULL, "alloc peano seen");
        if (!seen) return;
        for (x = 0; x < side; x++) {
            for (y = 0; y < side; y++) {
                d = curve_peano_xy2d(x, y, digits);
                ASSERT_TRUE(d < N, "peano index in range");
                if (d < N) {
                    ASSERT_TRUE(!seen[d], "peano injective");
                    seen[d] = 1;
                }
                curve_peano_d2xy(d, digits, &rx, &ry);
                ASSERT_INT_EQ(rx, x, "peano inverse x");
                ASSERT_INT_EQ(ry, y, "peano inverse y");
            }
        }
        for (i = 0; i < N; i++) {
            curve_peano_d2xy(i, digits, &x, &y);
            if (i > 0) {
                uint32_t dx = x > px ? x - px : px - x;
                uint32_t dy = y > py ? y - py : py - y;
                ASSERT_TRUE(dx + dy == 1u, "peano edge-adjacent");
            }
            px = x;
            py = y;
        }
        free(seen);
    }
}

static void test_peano_clamp_wrap(void)
{
    uint32_t x, y;
    ASSERT_INT_EQ(curve_peano_xy2d(1, 2, 0), curve_peano_xy2d(1, 2, 1),
                  "digits=0 clamps to 1");
    ASSERT_INT_EQ(curve_peano_xy2d(4, 5, 11), curve_peano_xy2d(4, 5, 10),
                  "digits=11 clamps to 10");
    ASSERT_INT_EQ(curve_peano_xy2d(3, 0, 1), curve_peano_xy2d(0, 0, 1),
                  "x wraps modulo 3^digits");
    curve_peano_d2xy(9, 1, &x, &y);
    ASSERT_TRUE(x == 0 && y == 0, "d wraps modulo 9^digits");
}

static void test_peano_sampled_high(void)
{
    int digits;
    uint32_t seed = 1u;
    for (digits = 6; digits <= 10; digits++) {
        uint32_t i, side = 1, k;
        for (k = 0; k < (uint32_t)digits; k++) side *= 3u;
        for (i = 0; i < 256; i++) {
            uint32_t x = (seed = seed * 1664525u + 1013904223u) % side;
            uint32_t y = (seed = seed * 1664525u + 1013904223u) % side;
            uint32_t d = curve_peano_xy2d(x, y, digits);
            uint32_t rx, ry;
            curve_peano_d2xy(d, digits, &rx, &ry);
            ASSERT_INT_EQ(rx, x, "peano high x");
            ASSERT_INT_EQ(ry, y, "peano high y");
        }
    }
}

static void test_mixed_uniform_matches_mortonnd(void)
{
    unsigned bits[8];
    uint32_t coords[8], back[8];
    uint64_t mixed[4], uniform[4];
    size_t i, n, b;
    for (n = 1; n <= 8; n++) {
        for (b = 1; b <= 8; b++) {
            for (i = 0; i < n; i++) {
                bits[i] = (unsigned)b;
                coords[i] = (uint32_t)(i * 17u + 3u);
            }
            ASSERT_INT_EQ(curve_key_bits_mixed(n, bits), n * b, "mixed bits sum");
            ASSERT_INT_EQ(curve_mortonnd_encode_words_mixed(
                              coords, n, bits, mixed, 4),
                          CURVE_OK, "mixed encode");
            ASSERT_INT_EQ(curve_mortonnd_encode_words(
                              coords, n, (unsigned)b, uniform, 4),
                          CURVE_OK, "uniform encode");
            ASSERT_TRUE(memcmp(mixed, uniform,
                               curve_key_words_required(n, (unsigned)b) *
                               sizeof(uint64_t)) == 0,
                        "uniform mixed == mortonnd");
            ASSERT_INT_EQ(curve_mortonnd_decode_words_mixed(
                              mixed, curve_key_words_required_mixed(n, bits),
                              n, bits, back, 8),
                          CURVE_OK, "mixed decode");
            for (i = 0; i < n; i++) {
                uint32_t mask = (b >= 32) ? 0xFFFFFFFFu : ((1u << b) - 1u);
                ASSERT_INT_EQ(back[i], coords[i] & mask, "mixed roundtrip");
            }
        }
    }
}

static void test_mixed_anisotropic_2d(void)
{
    unsigned bits[2] = { 8u, 4u };
    uint32_t x, y, rx, ry, coords[2], back[2];
    uint64_t key[1];
    uint8_t *seen;
    ASSERT_INT_EQ(curve_key_bits_mixed(2, bits), 12, "8+4=12");
    ASSERT_INT_EQ(curve_key_words_required_mixed(2, bits), 1, "fits one word");
    seen = (uint8_t *)calloc(1u << 12, 1);
    ASSERT_TRUE(seen != NULL, "alloc mixed seen");
    if (!seen) return;
    for (x = 0; x < 256; x++) {
        for (y = 0; y < 16; y++) {
            uint32_t d;
            coords[0] = x;
            coords[1] = y;
            ASSERT_INT_EQ(curve_mortonnd_encode_words_mixed(
                              coords, 2, bits, key, 1),
                          CURVE_OK, "mixed 2d encode");
            d = (uint32_t)(key[0] >> (64u - 12u));
            ASSERT_TRUE(d < (1u << 12), "mixed key in 12 bits");
            ASSERT_TRUE(!seen[d], "mixed injective on rectangle");
            seen[d] = 1;
            ASSERT_INT_EQ(curve_mortonnd_decode_words_mixed(
                              key, 1, 2, bits, back, 2),
                          CURVE_OK, "mixed 2d decode");
            rx = back[0];
            ry = back[1];
            ASSERT_INT_EQ(rx, x, "mixed rx");
            ASSERT_INT_EQ(ry, y, "mixed ry");
        }
    }
    free(seen);
    coords[0] = 0xFFFFFFFFu;
    coords[1] = 0xFFFFFFFFu;
    ASSERT_INT_EQ(curve_mortonnd_encode_words_mixed(coords, 2, bits, key, 1),
                  CURVE_OK, "mask extra bits");
    ASSERT_INT_EQ(curve_mortonnd_decode_words_mixed(key, 1, 2, bits, back, 2),
                  CURVE_OK, "decode masked");
    ASSERT_INT_EQ(back[0], 0xFFu, "x masked to 8");
    ASSERT_INT_EQ(back[1], 0xFu, "y masked to 4");
}

static void test_mixed_errors(void)
{
    unsigned bits[2] = { 8u, 4u };
    unsigned bad0[2] = { 0u, 4u };
    unsigned bad33[2] = { 8u, 33u };
    uint32_t coords[2] = { 1u, 2u };
    uint64_t key[1];
    uint32_t back[2];
    ASSERT_INT_EQ(curve_key_bits_mixed(0, bits), 0, "zero dims");
    ASSERT_INT_EQ(curve_key_bits_mixed(2, NULL), 0, "null bits");
    ASSERT_INT_EQ(curve_mortonnd_encode_words_mixed(coords, 2, bad0, key, 1),
                  CURVE_EINVAL, "bits=0");
    ASSERT_INT_EQ(curve_mortonnd_encode_words_mixed(coords, 2, bad33, key, 1),
                  CURVE_EINVAL, "bits=33");
    ASSERT_INT_EQ(curve_mortonnd_encode_words_mixed(NULL, 2, bits, key, 1),
                  CURVE_EINVAL, "null coords");
    ASSERT_INT_EQ(curve_mortonnd_encode_words_mixed(coords, 2, bits, key, 0),
                  CURVE_ENOSPACE, "short key");
    ASSERT_INT_EQ(curve_mortonnd_decode_words_mixed(key, 1, 2, bits, back, 1),
                  CURVE_ENOSPACE, "short coords");
}

int main(void)
{
    RUN_TEST(test_gmorton2_matches_gray_then_morton);
    RUN_TEST(test_gmorton2_64_roundtrip);
    RUN_TEST(test_gmorton3_roundtrip);
    RUN_TEST(test_gmortonnd_wide);
    RUN_TEST(test_peano_order1);
    RUN_TEST(test_peano_bijection_adjacency);
    RUN_TEST(test_peano_clamp_wrap);
    RUN_TEST(test_peano_sampled_high);
    RUN_TEST(test_mixed_uniform_matches_mortonnd);
    RUN_TEST(test_mixed_anisotropic_2d);
    RUN_TEST(test_mixed_errors);
    TEST_SUMMARY();
}
