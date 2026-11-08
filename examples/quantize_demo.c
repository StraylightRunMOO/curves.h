/* Quantize a 256-D vector and emit fused Morton / Hilbert keys. */

#include "curves.h"

#include <stdio.h>
#include <stdint.h>

int main(void)
{
    enum { N = 256, B = 4 };
    double v[N], mn[N], mx[N];
    uint32_t q[N], scratch[N];
    uint64_t morton[16], hilbert[16], fused_m[16], fused_h[16];
    size_t i, n_lo = 0, n_hi = 0, n_nan = 0;
    int mismatch = 0;
    curve_status st;

    printf("curves.h %s  isa=%s\n", CURVES_VERSION, curve_simd_isa());
    printf("256-D x %d-bit nearest, range [0, 1]\n", B);

    for (i = 0; i < (size_t)N; i++) {
        mn[i] = 0.0;
        mx[i] = 1.0;
        v[i] = (double)(i % 17) / 16.0;
    }
    v[0] = -0.25;
    v[1] = 1.5;
    {
        volatile double z = 0.0;
        v[2] = z / z; /* NaN */
    }

    st = curve_quantize_f64(v, N, mn, mx, B, CURVE_QUANT_NEAREST,
                            q, N, &n_lo, &n_hi, &n_nan);
    if (st != CURVE_OK) {
        fprintf(stderr, "quantize failed: %d\n", (int)st);
        return 1;
    }
    printf("saturated  lo=%zu  hi=%zu  nan=%zu\n", n_lo, n_hi, n_nan);

    st = curve_mortonnd_encode_words(q, N, B, morton, 16);
    if (st != CURVE_OK) {
        fprintf(stderr, "morton encode failed: %d\n", (int)st);
        return 1;
    }
    st = curve_hilbertnd_encode_words(q, N, B, hilbert, 16, scratch, N);
    if (st != CURVE_OK) {
        fprintf(stderr, "hilbert encode failed: %d\n", (int)st);
        return 1;
    }

    st = curve_mortonnd_encode_f64(v, N, mn, mx, B, CURVE_QUANT_NEAREST,
                                   fused_m, 16, scratch, N, NULL, NULL, NULL);
    if (st != CURVE_OK) {
        fprintf(stderr, "morton fused failed: %d\n", (int)st);
        return 1;
    }
    st = curve_hilbertnd_encode_f64(v, N, mn, mx, B, CURVE_QUANT_NEAREST,
                                    fused_h, 16, scratch, N, NULL, NULL, NULL);
    if (st != CURVE_OK) {
        fprintf(stderr, "hilbert fused failed: %d\n", (int)st);
        return 1;
    }

    for (i = 0; i < 16; i++) {
        if (morton[i] != fused_m[i] || hilbert[i] != fused_h[i])
            mismatch = 1;
    }
    printf("morton  word0=%016llx\n", (unsigned long long)morton[0]);
    printf("hilbert word0=%016llx\n", (unsigned long long)hilbert[0]);
    printf("fused==unfused: %s\n", mismatch ? "FAIL" : "ok");
    return mismatch;
}
