/* curves.h benchmarks — throughput of curve encode/decode and traversal. */

#define _POSIX_C_SOURCE 199309L

#include "curves.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static double now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

typedef struct {
    char name[48];
    double ns_per;
    double per_sec;
    uint64_t checksum;
} bench_row;

static bench_row g_rows[96];
static int g_nrows;
static int g_csv;

static void bench_record(const char *name, size_t iters, double dt, uint64_t checksum)
{
    bench_row *r;
    if (g_nrows >= (int)(sizeof(g_rows) / sizeof(g_rows[0]))) return;
    r = &g_rows[g_nrows++];
    snprintf(r->name, sizeof(r->name), "%s", name);
    r->ns_per = dt * 1e9 / (double)iters;
    r->per_sec = (double)iters / dt;
    r->checksum = checksum;
}

static void bench_print_ops(const char *name, size_t iters, double dt, uint64_t checksum)
{
    bench_record(name, iters, dt, checksum);
    if (!g_csv) {
        printf("%-28s %8.2f Mops/s  (%8.2f ns/op)\n", name,
               (double)iters / dt / 1e6, dt * 1e9 / (double)iters);
    }
}

#define BENCH(name, iters, ...) do {                                      \
    double t0 = now_sec();                                                \
    for (size_t _i = 0; _i < (iters); ++_i) { __VA_ARGS__ }               \
    double dt = now_sec() - t0;                                           \
    bench_print_ops((name), (iters), dt, (uint64_t)sink);                 \
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

int main(int argc, char **argv)
{
    int a;
    for (a = 1; a < argc; a++) {
        if (strcmp(argv[a], "--csv") == 0) g_csv = 1;
        else if (strcmp(argv[a], "--help") == 0) {
            printf("usage: bench_main [--csv]\n");
            return 0;
        }
    }

    if (!g_csv)
        printf("curves.h %s  isa=%s\n", CURVES_VERSION, curve_simd_isa());
    if (!g_csv) {
#if defined(CURVE_HAVE_BMI2)
        printf("BMI2 (PDEP/PEXT) scalar fast paths: enabled\n");
#else
        printf("BMI2 (PDEP/PEXT) scalar fast paths: portable fallback\n");
#endif
    }

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
    BENCH("gmorton2 encode", N, {
        sink += curve_gmorton2((uint32_t)_i, (uint32_t)(_i >> 8));
    });
    BENCH("peano xy2d (d=5)", N, {
        sink += curve_peano_xy2d((uint32_t)_i, (uint32_t)(_i >> 5), 5);
    });

    {
        uint32_t wcoords[256], wscratch[256], wout[256];
        uint64_t k2[1], k3[1], k16[1], k256_4[16], k256_8[32];
        size_t wi;
        const size_t N256 = N / 100;
        for (wi = 0; wi < 256; ++wi)
            wcoords[wi] = (uint32_t)wi * 2654435761u;
        curve_mortonnd_encode_words(wcoords, 2, 16, k2, 1);
        curve_mortonnd_encode_words(wcoords, 3, 10, k3, 1);
        curve_mortonnd_encode_words(wcoords, 16, 4, k16, 1);
        curve_mortonnd_encode_words(wcoords, 256, 4, k256_4, 16);
        curve_mortonnd_encode_words(wcoords, 256, 8, k256_8, 32);
        BENCH("mortonnd enc n=2 b=16", N, {
            wcoords[0] = (uint32_t)_i;
            curve_mortonnd_encode_words(wcoords, 2, 16, k2, 1);
            sink += k2[0];
        });
        BENCH("mortonnd dec n=2 b=16", N, {
            k2[0] = (uint64_t)_i << 32;
            curve_mortonnd_decode_words(k2, 1, 2, 16, wout, 256);
            sink += wout[0] + wout[1];
        });
        BENCH("mortonnd enc n=3 b=10", N, {
            wcoords[0] = (uint32_t)_i;
            curve_mortonnd_encode_words(wcoords, 3, 10, k3, 1);
            sink += k3[0];
        });
        BENCH("hilbertnd enc n=2 b=16", N, {
            wcoords[0] = (uint32_t)_i;
            curve_hilbertnd_encode_words(wcoords, 2, 16, k2, 1,
                                         wscratch, 256);
            sink += k2[0];
        });
        BENCH("hilbertnd dec n=2 b=16", N, {
            k2[0] = (uint64_t)_i << 32;
            curve_hilbertnd_decode_words(k2, 1, 2, 16, wout, 256,
                                         wscratch, 256);
            sink += wout[0] + wout[1];
        });
        BENCH("mortonnd enc n=16 b=4", N / 10, {
            wcoords[0] = (uint32_t)_i;
            curve_mortonnd_encode_words(wcoords, 16, 4, k16, 1);
            sink += k16[0];
        });
        BENCH("hilbertnd enc n=16 b=4", N / 10, {
            wcoords[0] = (uint32_t)_i;
            curve_hilbertnd_encode_words(wcoords, 16, 4, k16, 1,
                                         wscratch, 256);
            sink += k16[0];
        });
        BENCH("mortonnd enc n=256 b=4", N256, {
            wcoords[0] = (uint32_t)_i;
            curve_mortonnd_encode_words(wcoords, 256, 4, k256_4, 16);
            sink += k256_4[0] + k256_4[15];
        });
        BENCH("mortonnd dec n=256 b=4", N256, {
            k256_4[0] = (uint64_t)_i << 32;
            curve_mortonnd_decode_words(k256_4, 16, 256, 4, wout, 256);
            sink += wout[0] + wout[255];
        });
        BENCH("hilbertnd enc n=256 b=4", N256, {
            wcoords[0] = (uint32_t)_i;
            curve_hilbertnd_encode_words(wcoords, 256, 4, k256_4, 16,
                                         wscratch, 256);
            sink += k256_4[0] + k256_4[15];
        });
        BENCH("hilbertnd dec n=256 b=4", N256, {
            k256_4[0] = (uint64_t)_i << 32;
            curve_hilbertnd_decode_words(k256_4, 16, 256, 4, wout, 256,
                                         wscratch, 256);
            sink += wout[0] + wout[255];
        });
        BENCH("mortonnd enc n=256 b=8", N256, {
            wcoords[0] = (uint32_t)_i;
            curve_mortonnd_encode_words(wcoords, 256, 8, k256_8, 32);
            sink += k256_8[0] + k256_8[31];
        });
        BENCH("hilbertnd enc n=256 b=8", N256, {
            wcoords[0] = (uint32_t)_i;
            curve_hilbertnd_encode_words(wcoords, 256, 8, k256_8, 32,
                                         wscratch, 256);
            sink += k256_8[0] + k256_8[31];
        });
        {
            uint64_t ka[16], kb[16];
            size_t cp = 0;
            int m = 0;
            for (wi = 0; wi < 16; ++wi) {
                ka[wi] = 0x0123456789ABCDEFull + (uint64_t)wi;
                kb[wi] = ka[wi];
            }
            kb[15] ^= 1ull;
            BENCH("key_compare 16 words", N, {
                sink += (uint64_t)(curve_key_compare(ka, kb, 16) + 1);
            });
            BENCH("key_common_prefix 16w", N, {
                curve_key_common_prefix_bits(ka, kb, 1024, &cp);
                sink += cp;
            });
            BENCH("key_prefix_matches 16w", N, {
                curve_key_prefix_matches(ka, kb, 1024, 256, &m);
                sink += (uint64_t)m;
            });
        }
        {
            uint64_t bscratch[64], blo[64], bhi[64];
            int bunb[64];
            uint32_t baxes[16], bmin[2] = { 40u, 80u }, bmax[2] = { 90u, 140u };
            size_t nc = 0;
            BENCH("morton box 2D b=8", N / 100, {
                curve_morton_box_intervals(2, 8, bmin, bmax, blo, bhi, bunb,
                                           64, &nc, CURVE_RANGE_EXACT,
                                           bscratch, 64, baxes, 16);
                sink += nc;
            });
            BENCH("hilbert box 2D b=8", N / 100, {
                curve_hilbert_box_intervals(2, 8, bmin, bmax, blo, bhi, bunb,
                                            64, &nc, CURVE_RANGE_EXACT,
                                            bscratch, 64, baxes, 16);
                sink += nc;
            });
        }
        {
            double qv[256], qmin[256], qmax[256];
            uint32_t qout[256], qscratch[256];
            uint64_t qkey[16], qkey2[16];
            size_t qi;
            const size_t Nq = N / 100;
            for (qi = 0; qi < 256; ++qi) {
                qv[qi] = (double)qi * 0.01;
                qmin[qi] = 0.0;
                qmax[qi] = 4.0;
            }
            BENCH("quantize f64 n=256 b=4", Nq, {
                qv[0] = (double)(_i & 255u) * 0.01;
                curve_quantize_f64(qv, 256, qmin, qmax, 4, CURVE_QUANT_NEAREST,
                                   qout, 256, NULL, NULL, NULL);
                sink += qout[0] + qout[255];
            });
            BENCH("morton fused n=256 b=4", Nq, {
                qv[0] = (double)(_i & 255u) * 0.01;
                curve_mortonnd_encode_f64(qv, 256, qmin, qmax, 4,
                                          CURVE_QUANT_NEAREST,
                                          qkey, 16, qscratch, 256,
                                          NULL, NULL, NULL);
                sink += qkey[0] + qkey[15];
            });
            BENCH("morton unfused n=256 b=4", Nq, {
                qv[0] = (double)(_i & 255u) * 0.01;
                curve_quantize_f64(qv, 256, qmin, qmax, 4, CURVE_QUANT_NEAREST,
                                   qout, 256, NULL, NULL, NULL);
                curve_mortonnd_encode_words(qout, 256, 4, qkey2, 16);
                sink += qkey2[0] + qkey2[15];
            });
            BENCH("hilbert fused n=256 b=4", Nq, {
                qv[0] = (double)(_i & 255u) * 0.01;
                curve_hilbertnd_encode_f64(qv, 256, qmin, qmax, 4,
                                           CURVE_QUANT_NEAREST,
                                           qkey, 16, qscratch, 256,
                                           NULL, NULL, NULL);
                sink += qkey[0] + qkey[15];
            });
            BENCH("hilbert unfused n=256 b=4", Nq, {
                qv[0] = (double)(_i & 255u) * 0.01;
                curve_quantize_f64(qv, 256, qmin, qmax, 4, CURVE_QUANT_NEAREST,
                                   qout, 256, NULL, NULL, NULL);
                curve_hilbertnd_encode_words(qout, 256, 4, qkey2, 16,
                                             qscratch, 256);
                sink += qkey2[0] + qkey2[15];
            });
        }
    }

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
        bench_print_ops("morton2 batch", BN * repeats, dt, (uint64_t)sink);

        t0 = now_sec();
        for (size_t r = 0; r < repeats; ++r) {
            curve_morton2_decode_batch(bout, bx, by, BN);
            sink += bx[0] + by[0];
        }
        dt = now_sec() - t0;
        bench_print_ops("morton2 decode batch", BN * repeats, dt, (uint64_t)sink);

        t0 = now_sec();
        for (size_t r = 0; r < repeats; ++r) {
            curve_morton3_batch(bx, by, bz, bout64, BN);
            sink += bout64[0];
        }
        dt = now_sec() - t0;
        bench_print_ops("morton3 batch", BN * repeats, dt, (uint64_t)sink);

        t0 = now_sec();
        for (size_t r = 0; r < repeats; ++r) {
            size_t j;
            for (j = 0; j < BN; ++j)
                bout[j] = curve_hilbert_xy2d(bx[j], by[j], 16);
            sink += bout[0];
        }
        dt = now_sec() - t0;
        bench_print_ops("hilbert2 xy2d scalar", BN * repeats, dt, (uint64_t)sink);

        t0 = now_sec();
        for (size_t r = 0; r < repeats; ++r) {
            curve_hilbert_xy2d_batch(bx, by, 16, bout, BN);
            sink += bout[0];
        }
        dt = now_sec() - t0;
        bench_print_ops("hilbert2 xy2d batch", BN * repeats, dt, (uint64_t)sink);

        t0 = now_sec();
        for (size_t r = 0; r < repeats; ++r) {
            size_t j;
            for (j = 0; j < BN; ++j)
                curve_hilbert_d2xy(bout[j], 16, bx + j, by + j);
            sink += bx[0] + by[0];
        }
        dt = now_sec() - t0;
        bench_print_ops("hilbert2 d2xy scalar", BN * repeats, dt, (uint64_t)sink);

        t0 = now_sec();
        for (size_t r = 0; r < repeats; ++r) {
            curve_hilbert_d2xy_batch(bout, 16, bx, by, BN);
            sink += bx[0] + by[0];
        }
        dt = now_sec() - t0;
        bench_print_ops("hilbert2 d2xy batch", BN * repeats, dt, (uint64_t)sink);

        t0 = now_sec();
        for (size_t r = 0; r < repeats; ++r) {
            size_t j;
            for (j = 0; j < BN; ++j)
                bout64[j] = curve_hilbert2_64_xy2d(bx[j], by[j], 32);
            sink += bout64[0];
        }
        dt = now_sec() - t0;
        bench_print_ops("hilbert2_64 xy2d scalar", BN * repeats, dt, (uint64_t)sink);

        t0 = now_sec();
        for (size_t r = 0; r < repeats; ++r) {
            curve_hilbert2_64_xy2d_batch(bx, by, 32, bout64, BN);
            sink += bout64[0];
        }
        dt = now_sec() - t0;
        bench_print_ops("hilbert2_64 xy2d batch", BN * repeats, dt, (uint64_t)sink);

        t0 = now_sec();
        for (size_t r = 0; r < repeats; ++r) {
            curve_hilbert2_64_d2xy_batch(bout64, 32, bx, by, BN);
            sink += bx[0] + by[0];
        }
        dt = now_sec() - t0;
        bench_print_ops("hilbert2_64 d2xy batch", BN * repeats, dt, (uint64_t)sink);

        {
            curve_dispatch ds, dauto;
            unsigned bits2[2] = { 8u, 4u };
            uint32_t mc[2];
            uint64_t mk[1];
            curve_dispatch_init(&ds, CURVE_ISA_SCALAR);
            curve_dispatch_init(&dauto, 0);
            t0 = now_sec();
            for (size_t r = 0; r < repeats; ++r) {
                ds.hilbert_xy2d_batch(bx, by, 16, bout, BN);
                sink += bout[0];
            }
            dt = now_sec() - t0;
            bench_print_ops("hilbert2 xy2d disp.scalar", BN * repeats, dt,
                            (uint64_t)sink);
            t0 = now_sec();
            for (size_t r = 0; r < repeats; ++r) {
                dauto.hilbert_xy2d_batch(bx, by, 16, bout, BN);
                sink += bout[0];
            }
            dt = now_sec() - t0;
            bench_print_ops("hilbert2 xy2d disp.auto", BN * repeats, dt,
                            (uint64_t)sink);
            t0 = now_sec();
            for (size_t r = 0; r < repeats; ++r) {
                size_t j;
                for (j = 0; j < BN; ++j) {
                    mc[0] = bx[j];
                    mc[1] = by[j] & 15u;
                    curve_mortonnd_encode_words_mixed(mc, 2, bits2, mk, 1);
                    sink += mk[0];
                }
            }
            dt = now_sec() - t0;
            bench_print_ops("mixed morton 8+4", BN * repeats, dt, (uint64_t)sink);
        }
    }

    {
        double t0 = now_sec();
        for (int r = 0; r < 10; ++r)
            curve_gilbert2(1920, 1080, bench_count_cb2, &count);
        double dt = now_sec() - t0;
        bench_record("gilbert2 1920x1080", (size_t)count, dt, count);
        if (!g_csv)
            printf("%-28s %8.2f Mcells/s\n", "gilbert2 1920x1080",
                   (double)count / dt / 1e6);
    }
    {
        count = 0;
        double t0 = now_sec();
        for (int r = 0; r < 10; ++r)
            curve_gilbert3(64, 64, 64, bench_count_cb3, &count);
        double dt = now_sec() - t0;
        bench_record("gilbert3 64x64x64", (size_t)count, dt, count);
        if (!g_csv)
            printf("%-28s %8.2f Mcells/s\n", "gilbert3 64x64x64",
                   (double)count / dt / 1e6);
    }

    if (g_csv) {
        int i;
        printf("operation,median_ns,keys_per_sec,isa,checksum\n");
        for (i = 0; i < g_nrows; i++)
            printf("%s,%.4f,%.4f,%s,%llu\n",
                   g_rows[i].name, g_rows[i].ns_per, g_rows[i].per_sec,
                   curve_simd_isa(),
                   (unsigned long long)g_rows[i].checksum);
    }

    free(bx); free(by); free(bz); free(bout); free(bout64);
    return (int)(sink & 1) * 0;
}
