/* render_pgm — render space-filling curves to grayscale PGM images.
 *
 * Usage: render_pgm [outdir]
 * Writes: morton.pgm hilbert.pgm moore.pgm sierpinski.pgm gilbert.pgm
 * Pixel brightness = position along the curve (early = dark, late = light).
 */

#include "curves.h"

#include <stdio.h>
#include <stdlib.h>

static void write_pgm(const char *path, const uint8_t *img, int w, int h)
{
    FILE *f = fopen(path, "wb");
    if (!f) { perror(path); exit(1); }
    fprintf(f, "P5\n%d %d\n255\n", w, h);
    fwrite(img, 1, (size_t)(w * h), f);
    fclose(f);
    printf("wrote %s (%dx%d)\n", path, w, h);
}

static uint8_t *g_img;
static int g_w, g_h;
static uint64_t g_total;
static uint64_t g_i;

static void gilbert_paint(int x, int y, void *user)
{
    (void)user;
    g_img[y * g_w + x] = (uint8_t)(g_i * 255 / (g_total - 1));
    g_i++;
}

int main(int argc, char **argv)
{
    const char *dir = argc > 1 ? argv[1] : ".";
    char path[512];
    const int bits = 9;
    const int N = 1 << bits;                 /* 512x512 */
    const uint64_t total = (uint64_t)N * N;

    uint8_t *img = (uint8_t*)malloc((size_t)total);
    if (!img) return 1;

    for (uint32_t y = 0; y < (uint32_t)N; ++y)
        for (uint32_t x = 0; x < (uint32_t)N; ++x) {
            uint64_t i = (uint64_t)y * N + x;
            img[i] = (uint8_t)(curve_morton2(x, y) * 255 / (total - 1));
        }
    snprintf(path, sizeof path, "%s/morton.pgm", dir);
    write_pgm(path, img, N, N);

    for (uint32_t y = 0; y < (uint32_t)N; ++y)
        for (uint32_t x = 0; x < (uint32_t)N; ++x) {
            uint64_t i = (uint64_t)y * N + x;
            img[i] = (uint8_t)(curve_hilbert_xy2d(x, y, bits) * 255 / (total - 1));
        }
    snprintf(path, sizeof path, "%s/hilbert.pgm", dir);
    write_pgm(path, img, N, N);

    for (uint32_t y = 0; y < (uint32_t)N; ++y)
        for (uint32_t x = 0; x < (uint32_t)N; ++x) {
            uint64_t i = (uint64_t)y * N + x;
            img[i] = (uint8_t)(curve_moore_xy2d(x, y, bits) * 255 / (total - 1));
        }
    snprintf(path, sizeof path, "%s/moore.pgm", dir);
    write_pgm(path, img, N, N);

    /* Sierpiński key is non-dense in [0, 2*total); normalize by 2*total-1. */
    for (uint32_t y = 0; y < (uint32_t)N; ++y)
        for (uint32_t x = 0; x < (uint32_t)N; ++x) {
            uint64_t i = (uint64_t)y * N + x;
            img[i] = (uint8_t)(curve_sierpinski_xy2d(x, y, bits) * 255 /
                               (2 * total - 1));
        }
    snprintf(path, sizeof path, "%s/sierpinski.pgm", dir);
    write_pgm(path, img, N, N);

    free(img);

    /* Gilbert fills arbitrary rectangles, here 500x300 (not a power of two). */
    g_w = 500; g_h = 300;
    g_total = (uint64_t)g_w * g_h;
    g_i = 0;
    g_img = (uint8_t*)malloc((size_t)g_total);
    if (!g_img) return 1;
    curve_gilbert2(g_w, g_h, gilbert_paint, NULL);
    snprintf(path, sizeof path, "%s/gilbert.pgm", dir);
    write_pgm(path, g_img, g_w, g_h);
    free(g_img);

    return 0;
}
