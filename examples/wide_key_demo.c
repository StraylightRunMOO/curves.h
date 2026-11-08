/* Encode a 256-D vector at 4 bits/axis as 1024-bit Morton and Hilbert keys. */

#include "curves.h"

#include <stdio.h>
#include <stdint.h>

int main(void)
{
    uint32_t coords[256], scratch[256], back[256];
    uint64_t morton[16], hilbert[16];
    uint8_t bytes[128];
    size_t i;
    curve_status st;
    int mismatch = 0;

    for (i = 0; i < 256; i++)
        coords[i] = (uint32_t)(i & 15u); /* already quantized to 4 bits */

    printf("curves.h %s  isa=%s\n", CURVES_VERSION, curve_simd_isa());
    printf("256 dims x 4 bits = %zu key bits (%zu words, %zu bytes)\n",
           curve_key_bits(256, 4),
           curve_key_words_required(256, 4),
           curve_key_bytes_required(256, 4));

    st = curve_mortonnd_encode_words(coords, 256, 4, morton, 16);
    if (st != CURVE_OK) {
        fprintf(stderr, "morton encode failed: %d\n", (int)st);
        return 1;
    }
    st = curve_hilbertnd_encode_words(coords, 256, 4, hilbert, 16,
                                      scratch, 256);
    if (st != CURVE_OK) {
        fprintf(stderr, "hilbert encode failed: %d\n", (int)st);
        return 1;
    }

    printf("morton vs hilbert compare: %d\n",
           curve_key_compare(morton, hilbert, 16));

    st = curve_key_words_to_bytes(hilbert, 1024, bytes, 128);
    if (st != CURVE_OK) {
        fprintf(stderr, "bytes conversion failed: %d\n", (int)st);
        return 1;
    }
    printf("hilbert key (first 16 of 128 bytes):");
    for (i = 0; i < 16; i++)
        printf(" %02X", bytes[i]);
    printf(" ...\n");

    st = curve_hilbertnd_decode_words(hilbert, 16, 256, 4, back, 256,
                                      scratch, 256);
    if (st != CURVE_OK) {
        fprintf(stderr, "hilbert decode failed: %d\n", (int)st);
        return 1;
    }
    for (i = 0; i < 256; i++)
        if (back[i] != coords[i]) mismatch = 1;
    printf("hilbert roundtrip: %s\n", mismatch ? "FAIL" : "ok");
    return mismatch;
}
