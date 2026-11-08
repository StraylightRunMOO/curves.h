/* Decompose a 2D AABB into Morton and Hilbert key intervals. */

#include "curves.h"

#include <stdio.h>
#include <stdint.h>

int main(void)
{
    enum { BITS = 8, SIDE = 256 };
    uint32_t bmin[2] = { 40, 80 }, bmax[2] = { 90, 140 };
    uint64_t scratch[128], lo[256], hi[256];
    int unb[256];
    uint32_t axes[16];
    size_t n_m = 0, n_h = 0, i, need;
    curve_status st;

    printf("curves.h %s  isa=%s\n", CURVES_VERSION, curve_simd_isa());
    printf("grid %dx%d, box x[%u,%u] y[%u,%u]\n",
           SIDE, SIDE, bmin[0], bmax[0], bmin[1], bmax[1]);

    need = curve_box_intervals_scratch_words(2, BITS);
    if (need > 128u) {
        fprintf(stderr, "scratch too small: %zu\n", need);
        return 1;
    }

    st = curve_morton_box_intervals(2, BITS, bmin, bmax, lo, hi, unb, 256,
                                    &n_m, CURVE_RANGE_EXACT,
                                    scratch, 128, axes, 16);
    if (st != CURVE_OK) {
        fprintf(stderr, "morton box failed: %d\n", (int)st);
        return 1;
    }
    st = curve_hilbert_box_intervals(2, BITS, bmin, bmax, lo, hi, unb, 256,
                                     &n_h, CURVE_RANGE_EXACT,
                                     scratch, 128, axes, 16);
    if (st != CURVE_OK) {
        fprintf(stderr, "hilbert box failed: %d\n", (int)st);
        return 1;
    }

    printf("Morton intervals:  %zu\n", n_m);
    printf("Hilbert intervals: %zu\n", n_h);

    st = curve_hilbert_box_intervals(2, BITS, bmin, bmax, lo, hi, unb, 4,
                                     &n_h, CURVE_RANGE_OVERAPPROX,
                                     scratch, 128, axes, 16);
    if (st != CURVE_OK) {
        fprintf(stderr, "overapprox failed: %d\n", (int)st);
        return 1;
    }
    printf("Hilbert overapprox (cap=4): %zu spans\n", n_h);
    for (i = 0; i < n_h; i++)
        printf("  [%016llx, %s)\n",
               (unsigned long long)lo[i],
               unb[i] ? "+inf" : "next");
    return 0;
}
