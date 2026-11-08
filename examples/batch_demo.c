/* SIMD Hilbert 2D batch vs scalar, including a SIMD-body + tail length. */

#include "curves.h"

#include <stdio.h>
#include <stdint.h>

int main(void)
{
    enum { N = 17 };
    uint32_t x[N], y[N], d[N], rx[N], ry[N];
    uint64_t d64[N];
    size_t i;
    int mismatch = 0;

    printf("curves.h %s  compiled_isa=%s  simd_isa=%s\n",
           CURVES_VERSION, curve_compiled_isa(), curve_simd_isa());

    for (i = 0; i < (size_t)N; i++) {
        x[i] = (uint32_t)(i * 19u + 3u);
        y[i] = (uint32_t)(i * 11u + 7u);
    }

    curve_hilbert_xy2d_batch(x, y, 8, d, N);
    curve_hilbert_d2xy_batch(d, 8, rx, ry, N);
    curve_hilbert2_64_xy2d_batch(x, y, 20, d64, N);

    for (i = 0; i < (size_t)N; i++) {
        uint32_t sx, sy;
        if (d[i] != curve_hilbert_xy2d(x[i], y[i], 8)) mismatch = 1;
        curve_hilbert_d2xy(d[i], 8, &sx, &sy);
        if (rx[i] != sx || ry[i] != sy) mismatch = 1;
        if (d64[i] != curve_hilbert2_64_xy2d(x[i], y[i], 20)) mismatch = 1;
    }

    printf("n=%d (SIMD body + tail)  d[0]=%u  d[16]=%u\n",
           N, d[0], d[16]);
    printf("batch==scalar: %s\n", mismatch ? "FAIL" : "ok");
    return mismatch;
}
