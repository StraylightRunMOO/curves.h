/* curves.h benchmarks — throughput of curve encode/decode and traversal. */

#define _POSIX_C_SOURCE 199309L

#include "curves.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static double now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

#define BENCH(name, iters, ...) do {                                      \
    double t0 = now_sec();                                                \
    for (size_t _i = 0; _i < (iters); ++_i) { __VA_ARGS__ }               \
    double dt = now_sec() - t0;                                           \
    printf("%-28s %8.2f Mops/s  (%8.2f ns/op)\n", name,                   \
           (double)(iters) / dt / 1e6, dt * 1e9 / (double)(iters));       \
} while (0)

static void bench_count_cb2(int x, int y, void *user)
{
    (void)x; (void)y;
    (*(uint64_t*)user)++;
}

static void bench_count_cb3(int x, int y, int z, void *user)
{
    (void)x; (void)y; (void)z;
    (*(uint64_t*)user)++;
}

int main(void)
{
    printf("curves.h %s  isa=%s\n", CURVES_VERSION, curve_simd_isa());
#if defined(CURVE_HAVE_BMI2)
    printf("BMI2 (PDEP/PEXT) scalar fast paths: enabled\n");
#else
    printf("BMI2 (PDEP/PEXT) scalar fast paths: portable fallback\n");
#endif

    const size_t N = 1000000;
    volatile uint64_t sink = 0;
    uint64_t count = 0;

    enum { BN = 1 << 20 };
    uint32_t *bx = (uint32_t *)aligned_alloc(64, BN * sizeof(uint32_t));
    uint32_t *by = (uint32_t *)aligned_alloc(64, BN * sizeof(uint32_t));
    uint32_t *bz = (uint32_t *)aligned_alloc(64, BN * sizeof(uint32_t));
    uint32_t *bout = (uint32_t *)aligned_alloc(64, BN * sizeof(uint32_t));
    uint64_t *bout64 = (uint64_t *)aligned_alloc(64, BN * sizeof(uint64_t));
    if (!bx || !by || !bz || !bout || !bout64) {
        fprintf(stderr, "aligned_alloc failed\n");
        return 1;
    }
    for (size_t i = 0; i < BN; ++i) {
        bx[i] = (uint32_t)i * 2654435761u;
        by[i] = (uint32_t)i * 2246822519u + 13u;
        bz[i] = (uint32_t)i * 3266489917u + 29u;
    }

    BENCH("gray_encode", N * 10, { sink += curve_gray_encode(_i); });
    BENCH("gray_decode", N * 10, { sink += curve_gray_decode(_i); });

    BENCH("morton2 encode", N, {
        sink += curve_morton2((uint32_t)_i, (uint32_t)(_i >> 8));
    });
    BENCH("morton2 decode", N, {
        uint32_t x, y;
        curve_morton2_decode((uint32_t)_i * 1103515245u, &x, &y);
        sink += x + y;
    });
    BENCH("morton2_64 encode", N, {
        sink += curve_morton2_64((uint32_t)_i * 2654435761u,
                                 (uint32_t)(_i >> 3) * 2246822519u);
    });
    BENCH("morton2_64 decode", N, {
        uint32_t x, y;
        curve_morton2_64_decode((uint64_t)_i * 11400714819323198485ull, &x, &y);
        sink += x + y;
    });
    BENCH("morton3 encode", N, {
        sink += curve_morton3((uint32_t)_i & 1023u,
                              (uint32_t)(_i >> 10) & 1023u,
                              (uint32_t)(_i >> 20) & 1023u);
    });
    BENCH("morton3 decode", N, {
        uint32_t x, y, z;
        curve_morton3_decode((uint64_t)_i * 1099511628211ull, &x, &y, &z);
        sink += x + y + z;
    });

    BENCH("hilbert2 xy2d (16-bit)", N, {
        sink += curve_hilbert_xy2d((uint32_t)_i & 0xFFFFu,
                                   (uint32_t)(_i >> 16) & 0xFFFFu, 16);
    });
    BENCH("hilbert2 d2xy (16-bit)", N, {
        uint32_t x, y;
        curve_hilbert_d2xy((uint32_t)_i * 2654435761u, 16, &x, &y);
        sink += x + y;
    });
    BENCH("hilbert2_64 xy2d (32-bit)", N, {
        sink += curve_hilbert2_64_xy2d((uint32_t)_i * 2654435761u,
                                       (uint32_t)(_i >> 5) * 2246822519u, 32);
    });
    BENCH("hilbert2_64 d2xy (32-bit)", N, {
        uint32_t x, y;
        curve_hilbert2_64_d2xy((uint64_t)_i * 11400714819323198485ull, 32,
                               &x, &y);
        sink += x + y;
    });
    BENCH("hilbert3 xyz2d (10-bit)", N, {
        sink += curve_hilbert3_xyz2d((uint32_t)_i & 1023u,
                                     (uint32_t)(_i >> 10) & 1023u,
                                     (uint32_t)(_i >> 20) & 1023u, 10);
    });
    BENCH("hilbert3 d2xyz (10-bit)", N, {
        uint32_t x, y, z;
        curve_hilbert3_d2xyz((uint32_t)_i * 2654435761u, 10, &x, &y, &z);
        sink += x + y + z;
    });
    BENCH("hilbertnd xy2d (n=3,b=10)", N, {
        uint32_t c[3] = { (uint32_t)_i & 1023u, (uint32_t)(_i >> 10) & 1023u,
                          (uint32_t)(_i >> 20) & 1023u };
        sink += curve_hilbertnd_xy2d(c, 3, 10);
    });
    BENCH("hilbertnd d2xy (n=3,b=10)", N, {
        uint32_t c[3];
        curve_hilbertnd_d2xy((uint32_t)_i * 2654435761u, c, 3, 10);
        sink += c[0] + c[1] + c[2];
    });
    BENCH("moore xy2d (16-bit)", N, {
        sink += curve_moore_xy2d((uint32_t)_i & 0xFFFFu,
                                 (uint32_t)(_i >> 16) & 0xFFFFu, 16);
    });
    BENCH("moore d2xy (16-bit)", N, {
        uint32_t x, y;
        curve_moore_d2xy((uint32_t)_i * 2654435761u, 16, &x, &y);
        sink += x + y;
    });
    BENCH("sierpinski xy2d (16-bit)", N, {
        sink += curve_sierpinski_xy2d((uint32_t)_i & 0xFFFFu,
                                      (uint32_t)(_i >> 16) & 0xFFFFu, 16);
    });

    /* Batch: one call of BN elements, repeated to match ~N*something throughput.
     * Report Mops/s over BN * repeats elements. */
    {
        const size_t repeats = 16;
        double t0, dt;
        t0 = now_sec();
        for (size_t r = 0; r < repeats; ++r) {
            curve_morton2_batch(bx, by, bout, BN);
            sink += bout[0];
        }
        dt = now_sec() - t0;
        printf("%-28s %8.2f Mops/s  (%8.2f ns/op)\n", "morton2 batch",
               (double)(BN * repeats) / dt / 1e6,
               dt * 1e9 / (double)(BN * repeats));

        t0 = now_sec();
        for (size_t r = 0; r < repeats; ++r) {
            curve_morton2_decode_batch(bout, bx, by, BN);
            sink += bx[0] + by[0];
        }
        dt = now_sec() - t0;
        printf("%-28s %8.2f Mops/s  (%8.2f ns/op)\n", "morton2 decode batch",
               (double)(BN * repeats) / dt / 1e6,
               dt * 1e9 / (double)(BN * repeats));

        t0 = now_sec();
        for (size_t r = 0; r < repeats; ++r) {
            curve_morton3_batch(bx, by, bz, bout64, BN);
            sink += bout64[0];
        }
        dt = now_sec() - t0;
        printf("%-28s %8.2f Mops/s  (%8.2f ns/op)\n", "morton3 batch",
               (double)(BN * repeats) / dt / 1e6,
               dt * 1e9 / (double)(BN * repeats));

        t0 = now_sec();
        for (size_t r = 0; r < repeats; ++r) {
            curve_hilbert_xy2d_batch(bx, by, 16, bout, BN);
            sink += bout[0];
        }
        dt = now_sec() - t0;
        printf("%-28s %8.2f Mops/s  (%8.2f ns/op)\n", "hilbert2 xy2d batch",
               (double)(BN * repeats) / dt / 1e6,
               dt * 1e9 / (double)(BN * repeats));
    }

    {
        double t0 = now_sec();
        for (int r = 0; r < 10; ++r)
            curve_gilbert2(1920, 1080, bench_count_cb2, &count);
        double dt = now_sec() - t0;
        printf("%-28s %8.2f Mcells/s\n", "gilbert2 1920x1080",
               (double)count / dt / 1e6);
    }
    {
        count = 0;
        double t0 = now_sec();
        for (int r = 0; r < 10; ++r)
            curve_gilbert3(64, 64, 64, bench_count_cb3, &count);
        double dt = now_sec() - t0;
        printf("%-28s %8.2f Mcells/s\n", "gilbert3 64x64x64",
               (double)count / dt / 1e6);
    }

    free(bx); free(by); free(bz); free(bout); free(bout64);
    return (int)(sink & 1) * 0;
}
