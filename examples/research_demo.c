/* Gray-coded Morton, Peano 2D, mixed-bit keys, and caller-owned dispatch. */

#include "curves.h"

#include <stdio.h>
#include <stdint.h>

int main(void)
{
    uint32_t x, y, d, rx, ry;
    unsigned bits[2] = { 8u, 4u };
    uint32_t coords[2] = { 200u, 11u }, back[2];
    uint64_t key[1];
    curve_dispatch ds, dauto;
    uint32_t xs[8], ys[8], out_s[8], out_a[8];
    size_t i;
    int mismatch = 0;

    printf("curves.h %s  compiled=%s  detect=%s\n",
           CURVES_VERSION, curve_compiled_isa(),
           curve_isa_name(curve_isa_detect(0)));

    printf("gmorton2(3,5)=%u  morton2(3,5)=%u\n",
           curve_gmorton2(3, 5), curve_morton2(3, 5));
    curve_gmorton2_decode(curve_gmorton2(3, 5), &rx, &ry);
    printf("gmorton2 decode -> (%u,%u)\n", rx, ry);

    printf("peano order-1 path:");
    for (d = 0; d < 9; d++) {
        curve_peano_d2xy(d, 1, &x, &y);
        printf(" (%u,%u)", x, y);
    }
    printf("\n");
    printf("peano_xy2d(4,5,2)=%u\n", curve_peano_xy2d(4, 5, 2));

    curve_mortonnd_encode_words_mixed(coords, 2, bits, key, 1);
    curve_mortonnd_decode_words_mixed(key, 1, 2, bits, back, 2);
    printf("mixed 8+4 bits (%u,%u) key=0x%03x decode=(%u,%u)\n",
           coords[0], coords[1],
           (unsigned)(key[0] >> (64u - 12u)), back[0], back[1]);

    for (i = 0; i < 8; i++) {
        xs[i] = (uint32_t)(i * 19u + 3u);
        ys[i] = (uint32_t)(i * 11u + 7u);
    }
    curve_dispatch_init(&ds, CURVE_ISA_SCALAR);
    curve_dispatch_init(&dauto, 0);
    ds.morton2_batch(xs, ys, out_s, 8);
    dauto.morton2_batch(xs, ys, out_a, 8);
    for (i = 0; i < 8; i++)
        if (out_s[i] != out_a[i] || out_a[i] != curve_morton2(xs[i], ys[i]))
            mismatch = 1;
    printf("dispatch scalar=%s auto=%s match=%s\n",
           curve_selected_isa(&ds), curve_selected_isa(&dauto),
           mismatch ? "FAIL" : "ok");
    return mismatch;
}
