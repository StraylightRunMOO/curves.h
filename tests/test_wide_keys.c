#include "curves.h"
#include "test_common.h"
#include "ref_curves.h"

#include <stdlib.h>
#include <string.h>

static uint32_t mask_bits(unsigned bits)
{
    return bits >= 32u ? 0xFFFFFFFFu : ((1u << bits) - 1u);
}

static void fill_sentinel_u64(uint64_t *p, size_t n, uint64_t v)
{
    size_t i;
    for (i = 0; i < n; i++) p[i] = v;
}

static void fill_sentinel_u32(uint32_t *p, size_t n, uint32_t v)
{
    size_t i;
    for (i = 0; i < n; i++) p[i] = v;
}

static int words_equal(const uint64_t *a, const uint64_t *b, size_t n)
{
    size_t i;
    for (i = 0; i < n; i++)
        if (a[i] != b[i]) return 0;
    return 1;
}

static uint64_t left_align(uint64_t h, int K)
{
    if (K >= 64) return h;
    return h << (64 - K);
}

static void test_key_size_helpers(void)
{
    ASSERT_INT_EQ(curve_key_bits(256, 4), 1024, "256x4 bits");
    ASSERT_INT_EQ(curve_key_words_required(256, 4), 16, "256x4 words");
    ASSERT_INT_EQ(curve_key_bytes_required(256, 4), 128, "256x4 bytes");
    ASSERT_INT_EQ(curve_key_bits(256, 32), 8192, "256x32 bits");
    ASSERT_INT_EQ(curve_key_words_required(256, 32), 128, "256x32 words");
    ASSERT_INT_EQ(curve_key_bytes_required(256, 32), 1024, "256x32 bytes");
    ASSERT_INT_EQ(curve_key_bits(1, 1), 1, "1x1 bits");
    ASSERT_INT_EQ(curve_key_words_required(1, 1), 1, "1x1 words");
    ASSERT_INT_EQ(curve_key_bytes_required(1, 1), 1, "1x1 bytes");
    ASSERT_INT_EQ(curve_key_bits(3, 10), 30, "3x10 bits");
    ASSERT_INT_EQ(curve_key_words_required(2, 32), 1, "2x32 fits one word");
    ASSERT_INT_EQ(curve_key_bits(0, 4), 0, "zero dims invalid");
    ASSERT_INT_EQ(curve_key_bits(2, 0), 0, "zero bits invalid");
    ASSERT_INT_EQ(curve_key_bits(2, 33), 0, "bits>32 invalid");
    ASSERT_INT_EQ(curve_key_bits((size_t)CURVE_MAX_DIMS + 1u, 4), 0,
                  "too many dims invalid");
}

static void test_key_word_byte_roundtrip(void)
{
    uint64_t words[4], back[4];
    uint8_t bytes[32];
    int ok;

    /* K=16: one word, two bytes, known pattern. */
    words[0] = 0xABCD000000000000ull;
    ASSERT_INT_EQ(curve_key_words_to_bytes(words, 16, bytes, 2), CURVE_OK,
                  "words_to_bytes 16");
    ASSERT_INT_EQ(bytes[0], 0xAB, "byte0 is MSB");
    ASSERT_INT_EQ(bytes[1], 0xCD, "byte1");
    fill_sentinel_u64(back, 4, 0xDEADBEEFDEADBEEFull);
    ASSERT_INT_EQ(curve_key_bytes_to_words(bytes, 16, back, 4), CURVE_OK,
                  "bytes_to_words 16");
    ASSERT_INT_EQ(back[0], words[0], "16-bit word roundtrip");

    /* K=100: 2 words, 13 bytes, leftover nibble. */
    words[0] = 0x0123456789ABCDEFull;
    words[1] = 0xFEDCBA9870000000ull; /* 36 significant bits, 28 unused */
    ASSERT_INT_EQ(curve_key_words_to_bytes(words, 100, bytes, 13), CURVE_OK,
                  "words_to_bytes 100");
    ASSERT_INT_EQ(bytes[0], 0x01, "K=100 byte0");
    ASSERT_INT_EQ(bytes[12] & 0x0F, 0, "last-byte unused nibble is 0");
    fill_sentinel_u64(back, 4, ~0ull);
    ASSERT_INT_EQ(curve_key_bytes_to_words(bytes, 100, back, 2), CURVE_OK,
                  "bytes_to_words 100");
    ASSERT_TRUE(back[0] == words[0] && back[1] == words[1],
                "100-bit word/byte roundtrip");

    /* Noncanonical last byte is rejected. */
    bytes[12] |= 0x01;
    ASSERT_INT_EQ(curve_key_bytes_to_words(bytes, 100, back, 2), CURVE_EINVAL,
                  "noncanonical tail byte");

    /* Capacity failures leave the output untouched. */
    fill_sentinel_u64(back, 4, 0x1111111111111111ull);
    ASSERT_INT_EQ(curve_key_bytes_to_words(bytes, 16, back, 0), CURVE_ENOSPACE,
                  "short word cap");
    ASSERT_INT_EQ(back[0], 0x1111111111111111ull, "ENOSPACE does not write");

    fill_sentinel_u64(words, 1, 0x8000000000000000ull);
    bytes[0] = 0x5A;
    ASSERT_INT_EQ(curve_key_words_to_bytes(words, 1, bytes, 0), CURVE_ENOSPACE,
                  "short byte cap");
    ASSERT_INT_EQ(bytes[0], 0x5A, "ENOSPACE does not write bytes");

    ASSERT_INT_EQ(curve_key_words_to_bytes(NULL, 16, bytes, 2), CURVE_EINVAL,
                  "NULL words");
    ASSERT_INT_EQ(curve_key_words_to_bytes(words, 0, bytes, 2), CURVE_EINVAL,
                  "zero key_bits");

    /* Independent oracle vs conversion for a packed 4x4 vector. */
    {
        uint32_t X[4] = { 1, 2, 3, 4 };
        uint64_t packed[1];
        uint8_t b2[2], b3[2];
        ref_pack_axes_words(X, 4, 4, packed, 1);
        curve_key_words_to_bytes(packed, 16, b2, 2);
        /* second packing via u64 left-align */
        packed[0] = ref_pack_axes_u64(X, 4, 4) << (64 - 16);
        curve_key_words_to_bytes(packed, 16, b3, 2);
        ok = b2[0] == b3[0] && b2[1] == b3[1];
        ASSERT_TRUE(ok, "byte view matches u64-left-align oracle");
    }
}

static void test_key_compare(void)
{
    uint64_t a[2] = { 1ull, 2ull };
    uint64_t b[2] = { 1ull, 3ull };
    uint64_t c[2] = { 2ull, 0ull };
    ASSERT_INT_EQ(curve_key_compare(a, a, 2), 0, "equal");
    ASSERT_TRUE(curve_key_equal(a, a, 2), "equal helper");
    ASSERT_INT_EQ(curve_key_compare(a, b, 2), -1, "a < b on second word");
    ASSERT_INT_EQ(curve_key_compare(b, a, 2), 1, "b > a");
    ASSERT_INT_EQ(curve_key_compare(a, c, 2), -1, "first word dominates");
    ASSERT_INT_EQ(curve_key_compare(a, b, 0), 0, "zero words");
    /* Integer order for K<=64: left-aligned compare matches integer compare. */
    {
        uint64_t wa = left_align(3, 10);
        uint64_t wb = left_align(5, 10);
        ASSERT_TRUE(curve_key_compare(&wa, &wb, 1) < 0,
                    "left-align preserves integer order");
    }
}

static void test_checked_transpose(void)
{
    uint32_t axes[2] = { 3, 5 };
    uint32_t orig[2] = { 3, 5 };
    uint32_t one[1] = { 0xABCDu };
    ASSERT_INT_EQ(curve_axes_to_transpose(NULL, 2, 4), CURVE_EINVAL, "NULL");
    ASSERT_INT_EQ(curve_axes_to_transpose(axes, 0, 4), CURVE_EINVAL, "zero n");
    ASSERT_INT_EQ(curve_axes_to_transpose(axes, 2, 0), CURVE_EINVAL, "zero bits");
    ASSERT_INT_EQ(curve_axes_to_transpose(axes, 2, 33), CURVE_EINVAL, "bits 33");
    ASSERT_INT_EQ(curve_axes_to_transpose(axes, (size_t)CURVE_MAX_DIMS + 1u, 4),
                  CURVE_ERANGE, "too many dims");
    ASSERT_TRUE(axes[0] == orig[0] && axes[1] == orig[1],
                "failed transpose does not write");
    ASSERT_INT_EQ(curve_axes_to_transpose(one, 1, 16), CURVE_OK, "n=1 ok");
    ASSERT_INT_EQ(one[0], 0xABCDu, "n=1 is a no-op");

    /* n=2: transpose + pack matches legacy hilbertnd. */
    {
        uint32_t X[2] = { 3, 5 };
        uint64_t h, packed;
        int bits = 4, K = 8;
        curve_axes_to_transpose(X, 2, (unsigned)bits);
        packed = ref_pack_axes_u64(X, 2, bits);
        h = curve_hilbertnd_xy2d(orig, 2, bits);
        ASSERT_INT_EQ(packed, h, "transpose+pack matches hilbertnd");
        curve_transpose_to_axes(X, 2, (unsigned)bits);
        ASSERT_TRUE(X[0] == (orig[0] & mask_bits((unsigned)bits)) &&
                    X[1] == (orig[1] & mask_bits((unsigned)bits)),
                    "transpose inverse");
        (void)K;
    }
}

static void test_morton_failures(void)
{
    uint32_t coords[2] = { 1, 2 };
    uint64_t words[2];
    uint32_t out[2];
    fill_sentinel_u64(words, 2, 0xA5A5A5A5A5A5A5A5ull);
    ASSERT_INT_EQ(curve_mortonnd_encode_words(NULL, 2, 4, words, 2),
                  CURVE_EINVAL, "NULL coords");
    ASSERT_INT_EQ(words[0], 0xA5A5A5A5A5A5A5A5ull, "encode failure no write");
    ASSERT_INT_EQ(curve_mortonnd_encode_words(coords, 0, 4, words, 2),
                  CURVE_EINVAL, "zero dims");
    ASSERT_INT_EQ(curve_mortonnd_encode_words(coords, 2, 33, words, 2),
                  CURVE_EINVAL, "bits 33");
    ASSERT_INT_EQ(curve_mortonnd_encode_words(coords, 2, 4, words, 0),
                  CURVE_ENOSPACE, "short key cap");
    ASSERT_INT_EQ(curve_mortonnd_encode_words(coords,
                  (size_t)CURVE_MAX_DIMS + 1u, 4, words, 2),
                  CURVE_ERANGE, "too many dims");

    fill_sentinel_u32(out, 2, 0xDEADBEEFu);
    ASSERT_INT_EQ(curve_mortonnd_decode_words(NULL, 1, 2, 4, out, 2),
                  CURVE_EINVAL, "NULL key");
    ASSERT_INT_EQ(out[0], 0xDEADBEEFu, "decode failure no write");
}

static void test_morton_roundtrip_grid(void)
{
    static const size_t dims[] = { 1, 2, 3, 8, 16, 64, 65, 128, 256 };
    static const unsigned bitv[] = { 1, 2, 4, 8, 16, 32 };
    size_t di, bi, i;
    int ok = 1;
    uint32_t coords[256], back[256];
    uint64_t words[128];

    for (di = 0; di < sizeof(dims) / sizeof(dims[0]) && ok; ++di) {
        for (bi = 0; bi < sizeof(bitv) / sizeof(bitv[0]) && ok; ++bi) {
            size_t n = dims[di];
            unsigned bits = bitv[bi];
            size_t need = curve_key_words_required(n, bits);
            uint32_t m = mask_bits(bits);
            size_t sample, s;
            /* 256x32 is sampled; smaller configs get several vectors. */
            sample = (n >= 128 && bits >= 16) ? 8 : 32;
            for (s = 0; s < sample && ok; ++s) {
                for (i = 0; i < n; i++)
                    coords[i] = (uint32_t)(i * 2654435761u + s * 2246822519u);
                if (curve_mortonnd_encode_words(coords, n, bits, words, need)
                    != CURVE_OK) { ok = 0; break; }
                {
                    uint64_t refw[128];
                    uint32_t masked[256];
                    for (i = 0; i < n; i++) masked[i] = coords[i] & m;
                    ref_pack_axes_words(masked, (int)n, (int)bits, refw,
                                        (int)need);
                    if (!words_equal(words, refw, need)) { ok = 0; break; }
                }
                fill_sentinel_u32(back, n, 0xFFFFFFFFu);
                if (curve_mortonnd_decode_words(words, need, n, bits, back, n)
                    != CURVE_OK) { ok = 0; break; }
                for (i = 0; i < n; i++)
                    if (back[i] != (coords[i] & m)) { ok = 0; break; }
            }
        }
    }
    ASSERT_TRUE(ok, "mortonnd roundtrip vs ref pack, many (n,bits)");
}

static void test_morton_exhaustive_small(void)
{
    int bits, ok = 1;
    {
        unsigned x, y;
        for (bits = 1; bits <= 4 && ok; ++bits) {
            uint32_t N = 1u << bits;
            for (x = 0; x < N && ok; ++x)
                for (y = 0; y < N && ok; ++y) {
                    uint32_t c[2] = { x, y }, b[2];
                    uint64_t words[1], refw[1];
                    if (curve_mortonnd_encode_words(c, 2, (unsigned)bits,
                                                    words, 1) != CURVE_OK)
                    { ok = 0; break; }
                    ref_pack_axes_words(c, 2, bits, refw, 1);
                    if (words[0] != refw[0]) { ok = 0; break; }
                    if (curve_mortonnd_decode_words(words, 1, 2,
                            (unsigned)bits, b, 2) != CURVE_OK) { ok = 0; break; }
                    if (b[0] != x || b[1] != y) { ok = 0; break; }
                }
        }
    }
    ASSERT_TRUE(ok, "morton 2D exhaustive bits 1..4");

    ok = 1;
    {
        unsigned x, y, z;
        for (x = 0; x < 8 && ok; ++x)
            for (y = 0; y < 8 && ok; ++y)
                for (z = 0; z < 8 && ok; ++z) {
                    uint32_t c[3] = { x, y, z }, b[3];
                    uint64_t words[1], refw[1];
                    curve_mortonnd_encode_words(c, 3, 3, words, 1);
                    ref_pack_axes_words(c, 3, 3, refw, 1);
                    if (words[0] != refw[0]) { ok = 0; break; }
                    curve_mortonnd_decode_words(words, 1, 3, 3, b, 3);
                    if (b[0] != x || b[1] != y || b[2] != z) { ok = 0; break; }
                }
    }
    ASSERT_TRUE(ok, "morton 3D exhaustive bits=3");
}

static void test_morton_legacy_mismatch_is_documented(void)
{
    uint32_t c[2] = { 1, 0 };
    uint64_t words[1];
    curve_mortonnd_encode_words(c, 2, 1, words, 1);
    /* Wide packing: MSB-first, axis 0 first → 0b10 left-aligned.
     * Legacy morton2(1,0) == 1 (axis 0 on even LSB). */
    ASSERT_INT_EQ(curve_morton2(1, 0), 1, "legacy morton2(1,0)");
    ASSERT_TRUE(words[0] != 1ull, "wide morton n=2 is not curve_morton2");
    ASSERT_INT_EQ(words[0], 0x8000000000000000ull,
                  "wide morton(1,0) bits=1 is 0b10 left-aligned");
}

static void test_morton_256x4_known(void)
{
    uint32_t coords[256], back[256];
    uint64_t words[16], refw[16];
    uint8_t bytes[128], bytes2[128];
    size_t i;
    int ok;

    for (i = 0; i < 256; i++) coords[i] = 0;
    ASSERT_INT_EQ(curve_mortonnd_encode_words(coords, 256, 4, words, 16),
                  CURVE_OK, "256x4 zeros encode");
    for (i = 0; i < 16; i++)
        if (words[i] != 0) { ASSERT_TRUE(0, "zero vector → zero key"); return; }

    for (i = 0; i < 256; i++) coords[i] = 0xFu;
    curve_mortonnd_encode_words(coords, 256, 4, words, 16);
    ref_pack_axes_words(coords, 256, 4, refw, 16);
    ASSERT_TRUE(words_equal(words, refw, 16), "256x4 all-ones vs ref");
    curve_mortonnd_decode_words(words, 16, 256, 4, back, 256);
    ok = 1;
    for (i = 0; i < 256; i++)
        if (back[i] != 0xFu) ok = 0;
    ASSERT_TRUE(ok, "256x4 all-ones decode");

    for (i = 0; i < 256; i++) coords[i] = (uint32_t)(i & 15u);
    curve_mortonnd_encode_words(coords, 256, 4, words, 16);
    curve_key_words_to_bytes(words, 1024, bytes, 128);
    ref_pack_axes_words(coords, 256, 4, refw, 16);
    curve_key_words_to_bytes(refw, 1024, bytes2, 128);
    ok = 1;
    for (i = 0; i < 128; i++)
        if (bytes[i] != bytes2[i]) ok = 0;
    ASSERT_TRUE(ok, "256x4 byte serialization matches independent packer");

    /* walking-1 on axis 0 */
    for (i = 0; i < 256; i++) coords[i] = 0;
    coords[0] = 1;
    curve_mortonnd_encode_words(coords, 256, 4, words, 16);
    curve_mortonnd_decode_words(words, 16, 256, 4, back, 256);
    ok = (back[0] == 1);
    for (i = 1; i < 256; i++)
        if (back[i] != 0) ok = 0;
    ASSERT_TRUE(ok, "256x4 walking-1 axis0");

    /* 256x8 smoke */
    {
        uint64_t w32[32];
        uint32_t c8[256], b8[256];
        for (i = 0; i < 256; i++) c8[i] = (uint32_t)(i & 255u);
        ASSERT_INT_EQ(curve_mortonnd_encode_words(c8, 256, 8, w32, 32),
                      CURVE_OK, "256x8 encode");
        ASSERT_INT_EQ(curve_mortonnd_decode_words(w32, 32, 256, 8, b8, 256),
                      CURVE_OK, "256x8 decode");
        ok = 1;
        for (i = 0; i < 256; i++)
            if (b8[i] != c8[i]) ok = 0;
        ASSERT_TRUE(ok, "256x8 roundtrip");
    }
}

static void test_morton_noncanonical(void)
{
    uint32_t c[2] = { 1, 0 }, out[2];
    uint64_t words[1];
    curve_mortonnd_encode_words(c, 2, 1, words, 1);
    words[0] |= 1ull; /* pollute unused low bits */
    fill_sentinel_u32(out, 2, 9);
    ASSERT_INT_EQ(curve_mortonnd_decode_words(words, 1, 2, 1, out, 2),
                  CURVE_EINVAL, "noncanonical tail");
    ASSERT_INT_EQ(out[0], 9, "noncanonical decode does not write");
    ASSERT_INT_EQ(curve_mortonnd_decode_words(words, 2, 2, 1, out, 2),
                  CURVE_EINVAL, "wrong word count");
}

static void test_hilbert_failures(void)
{
    uint32_t coords[2] = { 1, 2 }, scratch[2], out[2];
    uint64_t words[2];
    fill_sentinel_u64(words, 2, 0xA5A5A5A5A5A5A5A5ull);
    ASSERT_INT_EQ(curve_hilbertnd_encode_words(coords, 2, 4, words, 2, NULL, 2),
                  CURVE_EINVAL, "NULL scratch");
    ASSERT_INT_EQ(words[0], 0xA5A5A5A5A5A5A5A5ull, "no write on NULL scratch");
    ASSERT_INT_EQ(curve_hilbertnd_encode_words(coords, 2, 4, words, 2,
                                              scratch, 1),
                  CURVE_ENOSPACE, "short scratch");
    ASSERT_INT_EQ(curve_hilbertnd_encode_words(NULL, 2, 4, words, 2,
                                              scratch, 2),
                  CURVE_EINVAL, "NULL coords");
    fill_sentinel_u32(out, 2, 7);
    ASSERT_INT_EQ(curve_hilbertnd_decode_words(words, 1, 2, 4, out, 2,
                                              NULL, 2),
                  CURVE_EINVAL, "NULL decode scratch");
    ASSERT_INT_EQ(out[0], 7, "decode failure no write");
}

static void test_hilbert_legacy_parity(void)
{
    int n, bits, ok = 1;
    for (n = 2; n <= 8 && ok; ++n) {
        for (bits = 1; bits <= 8 && ok; ++bits) {
            int K = n * bits;
            uint64_t total, step, d;
            if (K > 24) break;
            total = 1ull << K;
            step = total > 4096 ? total / 4096 : 1;
            for (d = 0; d < total && ok; d += step) {
                uint32_t coords[8], scratch[8], back[8];
                uint64_t words[1];
                uint64_t h;
                int i;
                for (i = 0; i < n; i++) coords[i] = 0;
                curve_hilbertnd_d2xy(d, coords, n, bits);
                h = curve_hilbertnd_xy2d(coords, n, bits);
                if (curve_hilbertnd_encode_words(coords, (size_t)n,
                        (unsigned)bits, words, 1, scratch, (size_t)n)
                    != CURVE_OK) { ok = 0; break; }
                if (words[0] != left_align(h, K)) { ok = 0; break; }
                if (n == 2 &&
                    words[0] != left_align(curve_hilbert_xy2d(coords[0],
                                                             coords[1], bits),
                                           K))
                { ok = 0; break; }
                fill_sentinel_u32(back, (size_t)n, 0xFFFFFFFFu);
                if (curve_hilbertnd_decode_words(words, 1, (size_t)n,
                        (unsigned)bits, back, (size_t)n, scratch, (size_t)n)
                    != CURVE_OK) { ok = 0; break; }
                for (i = 0; i < n; i++)
                    if (back[i] != coords[i]) { ok = 0; break; }
            }
        }
    }
    ASSERT_TRUE(ok, "wide hilbert matches legacy hilbertnd for n*bits<=24");
}

static void test_hilbert_exhaustive_adjacency(void)
{
    int bits, ok = 1, adj = 1, bij = 1;
    for (bits = 1; bits <= 5 && ok; ++bits) {
        uint32_t N = 1u << bits;
        uint64_t total = (uint64_t)N * N;
        uint8_t *seen = (uint8_t *)calloc((size_t)total, 1);
        uint32_t px = 0, py = 0;
        uint64_t d;
        if (!seen) { ASSERT_TRUE(0, "malloc"); return; }
        for (d = 0; d < total; ++d) {
            uint32_t c[2], scratch[2];
            uint64_t words[1];
            words[0] = left_align(d, 2 * bits);
            if (curve_hilbertnd_decode_words(words, 1, 2, (unsigned)bits,
                                             c, 2, scratch, 2) != CURVE_OK)
            { ok = 0; break; }
            if (c[0] >= N || c[1] >= N || seen[(uint64_t)c[1] * N + c[0]])
            { bij = 0; break; }
            seen[(uint64_t)c[1] * N + c[0]] = 1;
            if (d > 0) {
                uint32_t dx = c[0] > px ? c[0] - px : px - c[0];
                uint32_t dy = c[1] > py ? c[1] - py : py - c[1];
                if (dx + dy != 1) adj = 0;
            }
            px = c[0];
            py = c[1];
            {
                uint64_t enc[1];
                uint32_t sc[2];
                curve_hilbertnd_encode_words(c, 2, (unsigned)bits, enc, 1,
                                             sc, 2);
                if (enc[0] != words[0]) { ok = 0; break; }
            }
        }
        free(seen);
    }
    ASSERT_TRUE(ok && bij, "wide hilbert 2D bijective bits 1..5");
    ASSERT_TRUE(adj, "wide hilbert 2D edge-adjacent bits 1..5");

    /* 3D bits=1..3: Manhattan 1. */
    ok = 1;
    adj = 1;
    bij = 1;
    for (bits = 1; bits <= 3 && ok; ++bits) {
        uint32_t N = 1u << bits;
        uint64_t total = (uint64_t)N * N * N;
        uint8_t *seen = (uint8_t *)calloc((size_t)total, 1);
        uint32_t px = 0, py = 0, pz = 0;
        uint64_t d;
        if (!seen) { ASSERT_TRUE(0, "malloc"); return; }
        for (d = 0; d < total; ++d) {
            uint32_t c[3], scratch[3];
            uint64_t words[1];
            words[0] = left_align(d, 3 * bits);
            curve_hilbertnd_decode_words(words, 1, 3, (unsigned)bits, c, 3,
                                         scratch, 3);
            {
                uint64_t idx = ((uint64_t)c[2] * N + c[1]) * N + c[0];
                if (c[0] >= N || c[1] >= N || c[2] >= N || seen[idx])
                { bij = 0; break; }
                seen[idx] = 1;
            }
            if (d > 0) {
                uint32_t dx = c[0] > px ? c[0] - px : px - c[0];
                uint32_t dy = c[1] > py ? c[1] - py : py - c[1];
                uint32_t dz = c[2] > pz ? c[2] - pz : pz - c[2];
                if (dx + dy + dz != 1) adj = 0;
            }
            px = c[0]; py = c[1]; pz = c[2];
        }
        free(seen);
    }
    ASSERT_TRUE(ok && bij, "wide hilbert 3D bijective bits 1..3");
    ASSERT_TRUE(adj, "wide hilbert 3D face-adjacent bits 1..3");

    /* n=4 bits=2: 256 points. */
    ok = 1;
    bij = 1;
    adj = 1;
    {
        uint8_t seen[256];
        uint32_t prev[4];
        uint64_t d;
        memset(seen, 0, sizeof(seen));
        for (d = 0; d < 256; ++d) {
            uint32_t c[4], scratch[4];
            uint64_t words[1];
            int i, manh = 0;
            words[0] = left_align(d, 8);
            curve_hilbertnd_decode_words(words, 1, 4, 2, c, 4, scratch, 4);
            {
                unsigned code = 0;
                for (i = 0; i < 4; i++) {
                    if (c[i] > 3) { bij = 0; break; }
                    code = (code << 2) | c[i];
                }
                if (seen[code]) bij = 0;
                seen[code] = 1;
            }
            if (d > 0) {
                for (i = 0; i < 4; i++) {
                    uint32_t dx = c[i] > prev[i] ? c[i] - prev[i] : prev[i] - c[i];
                    manh += (int)dx;
                }
                if (manh != 1) adj = 0;
            }
            for (i = 0; i < 4; i++) prev[i] = c[i];
        }
    }
    ASSERT_TRUE(ok && bij, "wide hilbert n=4 bits=2 bijective");
    ASSERT_TRUE(adj, "wide hilbert n=4 bits=2 adjacent");
}

static void test_hilbert_256_and_alias(void)
{
    uint32_t coords[256], back[256];
    uint64_t words[16], w8[32];
    uint8_t bytes[128], bytes2[128];
    size_t i;
    int ok;

    for (i = 0; i < 256; i++) coords[i] = (uint32_t)(i & 15u);
    ASSERT_INT_EQ(curve_hilbertnd_encode_words(coords, 256, 4, words, 16,
                                               back, 256),
                  CURVE_OK, "256x4 hilbert encode");
    ASSERT_INT_EQ(curve_hilbertnd_decode_words(words, 16, 256, 4, back, 256,
                                               coords, 256),
                  CURVE_OK, "256x4 hilbert decode (scratch aliases coords)");
    /* decode overwrote coords via alias; recover from back vs original mask */
    {
        uint32_t orig[256];
        for (i = 0; i < 256; i++) orig[i] = (uint32_t)(i & 15u);
        ok = 1;
        for (i = 0; i < 256; i++)
            if (back[i] != orig[i]) ok = 0;
        ASSERT_TRUE(ok, "256x4 hilbert roundtrip");
    }

    for (i = 0; i < 256; i++) coords[i] = (uint32_t)(i & 15u);
    curve_hilbertnd_encode_words(coords, 256, 4, words, 16, back, 256);
    curve_key_words_to_bytes(words, 1024, bytes, 128);
    /* Independent: mask, Skilling via public transpose, ref pack. */
    {
        uint32_t X[256];
        uint64_t refw[16];
        for (i = 0; i < 256; i++) X[i] = coords[i] & 15u;
        curve_axes_to_transpose(X, 256, 4);
        ref_pack_axes_words(X, 256, 4, refw, 16);
        curve_key_words_to_bytes(refw, 1024, bytes2, 128);
        ok = 1;
        for (i = 0; i < 128; i++)
            if (bytes[i] != bytes2[i]) ok = 0;
        ASSERT_TRUE(ok, "256x4 hilbert bytes match transpose+ref pack");
        ASSERT_TRUE(words_equal(words, refw, 16), "256x4 hilbert words vs ref");
    }

    /* in-place: scratch == coords */
    for (i = 0; i < 256; i++) coords[i] = (uint32_t)((i * 3u) & 15u);
    {
        uint32_t saved[256];
        for (i = 0; i < 256; i++) saved[i] = coords[i];
        curve_hilbertnd_encode_words(coords, 256, 4, words, 16, coords, 256);
        curve_hilbertnd_decode_words(words, 16, 256, 4, coords, 256, coords, 256);
        ok = 1;
        for (i = 0; i < 256; i++)
            if (coords[i] != saved[i]) ok = 0;
        ASSERT_TRUE(ok, "in-place scratch==coords roundtrip");
    }

    for (i = 0; i < 256; i++) coords[i] = (uint32_t)(i & 255u);
    ASSERT_INT_EQ(curve_hilbertnd_encode_words(coords, 256, 8, w8, 32,
                                               back, 256),
                  CURVE_OK, "256x8 encode");
    ASSERT_INT_EQ(curve_hilbertnd_decode_words(w8, 32, 256, 8, back, 256,
                                               coords, 256),
                  CURVE_OK, "256x8 decode");
    ok = 1;
    for (i = 0; i < 256; i++)
        if (back[i] != (uint32_t)(i & 255u)) ok = 0;
    ASSERT_TRUE(ok, "256x8 hilbert roundtrip");
}

static void test_hilbert_n2_bits32(void)
{
    int i, ok = 1;
    for (i = 0; i < 20000 && ok; ++i) {
        uint32_t coords[2] = { (uint32_t)i * 2654435761u,
                               (uint32_t)i * 2246822519u + 7u };
        uint32_t scratch[2], back[2];
        uint64_t words[1];
        uint64_t h = curve_hilbertnd_xy2d(coords, 2, 32);
        curve_hilbertnd_encode_words(coords, 2, 32, words, 1, scratch, 2);
        if (words[0] != h) { ok = 0; break; }
        curve_hilbertnd_decode_words(words, 1, 2, 32, back, 2, scratch, 2);
        if (back[0] != coords[0] || back[1] != coords[1]) ok = 0;
    }
    ASSERT_TRUE(ok, "wide hilbert n=2 bits=32 matches legacy");
}

static void test_hilbert_n1_identity(void)
{
    uint32_t c[1] = { 0x12345678u }, b[1], scratch[1];
    uint64_t words[1];
    curve_hilbertnd_encode_words(c, 1, 32, words, 1, scratch, 1);
    ASSERT_INT_EQ(words[0], 0x1234567800000000ull,
                  "1D hilbert is the axis, left-aligned");
    curve_mortonnd_encode_words(c, 1, 32, words, 1);
    ASSERT_INT_EQ(words[0], 0x1234567800000000ull,
                  "1D morton is the axis, left-aligned");
    curve_hilbertnd_decode_words(words, 1, 1, 32, b, 1, scratch, 1);
    ASSERT_INT_EQ(b[0], 0x12345678u, "1D decode");
}

int main(void)
{
    RUN_TEST(test_key_size_helpers);
    RUN_TEST(test_key_word_byte_roundtrip);
    RUN_TEST(test_key_compare);
    RUN_TEST(test_checked_transpose);
    RUN_TEST(test_morton_failures);
    RUN_TEST(test_morton_roundtrip_grid);
    RUN_TEST(test_morton_exhaustive_small);
    RUN_TEST(test_morton_legacy_mismatch_is_documented);
    RUN_TEST(test_morton_256x4_known);
    RUN_TEST(test_morton_noncanonical);
    RUN_TEST(test_hilbert_failures);
    RUN_TEST(test_hilbert_legacy_parity);
    RUN_TEST(test_hilbert_exhaustive_adjacency);
    RUN_TEST(test_hilbert_256_and_alias);
    RUN_TEST(test_hilbert_n2_bits32);
    RUN_TEST(test_hilbert_n1_identity);
    TEST_SUMMARY();
}
