/* locality_demo — compare the locality of grid traversal orderings.
 *
 * For each ordering of a 2^bits x 2^bits grid, report:
 *   - mean step length between consecutive indices (1.0 = perfectly local)
 *   - mean bounding-box side of windows of 64 consecutive indices
 *     (smaller = better clustering; sqrt(64) = 8 is ideal for a square)
 */

#include "curves.h"

#include <stdio.h>
#include <stdlib.h>

typedef struct {
    const char *name;
    void (*decode)(uint64_t d, uint32_t *x, uint32_t *y, int bits);
} ordering;

static void dec_rowmajor(uint64_t d, uint32_t *x, uint32_t *y, int bits)
{
    *x = (uint32_t)(d & ((1u << bits) - 1u));
    *y = (uint32_t)(d >> bits);
}

static void dec_morton(uint64_t d, uint32_t *x, uint32_t *y, int bits)
{
    (void)bits;
    curve_morton2_decode((uint32_t)d, x, y);
}

static void dec_hilbert(uint64_t d, uint32_t *x, uint32_t *y, int bits)
{
    curve_hilbert_d2xy((uint32_t)d, bits, x, y);
}

static void dec_moore(uint64_t d, uint32_t *x, uint32_t *y, int bits)
{
    curve_moore_d2xy((uint32_t)d, bits, x, y);
}

int main(void)
{
    const int bits = 10;
    const uint32_t N = 1u << bits;
    const uint64_t total = (uint64_t)N * N;
    const uint64_t window = 64;

    const ordering orderings[] = {
        { "row-major", dec_rowmajor },
        { "morton (z-order)", dec_morton },
        { "hilbert", dec_hilbert },
        { "moore", dec_moore },
    };

    printf("grid: %ux%u (%llu cells), window = %llu\n\n", N, N,
           (unsigned long long)total, (unsigned long long)window);
    printf("%-18s %14s %16s\n", "ordering", "mean step", "mean bbox side");
    printf("%-18s %14s %16s\n", "--------", "---------", "--------------");

    for (size_t o = 0; o < sizeof(orderings) / sizeof(orderings[0]); ++o) {
        double step_sum = 0.0, bbox_sum = 0.0;
        uint64_t nbbox = 0;
        uint32_t px = 0, py = 0;
        uint32_t minx = 0, miny = 0, maxx = 0, maxy = 0;

        for (uint64_t d = 0; d < total; ++d) {
            uint32_t x, y;
            orderings[o].decode(d, &x, &y, bits);
            if (d == 0) {
                minx = maxx = x; miny = maxy = y;
            } else {
                double dx = (double)x - (double)px;
                double dy = (double)y - (double)py;
                step_sum += dx * dx + dy * dy;      /* squared step */
                if (x < minx) minx = x;
                if (x > maxx) maxx = x;
                if (y < miny) miny = y;
                if (y > maxy) maxy = y;
            }
            if ((d + 1) % window == 0) {
                uint64_t bw = maxx - minx + 1, bh = maxy - miny + 1;
                bbox_sum += (double)(bw > bh ? bw : bh);
                nbbox++;
                minx = maxx = x; miny = maxy = y;
            }
            px = x; py = y;
        }
        printf("%-18s %14.3f %16.2f\n", orderings[o].name,
               step_sum / (double)(total - 1), bbox_sum / (double)nbbox);
    }
    return 0;
}
