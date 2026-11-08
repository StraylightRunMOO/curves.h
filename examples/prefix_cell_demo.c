/* Prefix search and Hilbert cell navigation on a small 2D grid. */

#include "curves.h"

#include <stdio.h>
#include <stdlib.h>

static int cmp_u64(const void *a, const void *b)
{
    return curve_key_compare((const uint64_t *)a, (const uint64_t *)b, 1);
}

int main(void)
{
    enum { N = 16, BITS = 4 };
    uint32_t coords[2], scratch[2], mins[2], maxs[2];
    uint64_t keys[N * N], scratch_range[2], child_words[1];
    curve_cell cell, child;
    size_t i, x, y, begin, end, nchild;
    uint64_t prefix;
    curve_status st;

    printf("curves.h %s  isa=%s\n", CURVES_VERSION, curve_simd_isa());

    for (y = 0; y < N; y++) {
        for (x = 0; x < N; x++) {
            coords[0] = (uint32_t)x;
            coords[1] = (uint32_t)y;
            curve_hilbertnd_encode_words(coords, 2, BITS,
                                         &keys[y * N + x], 1,
                                         scratch, 2);
        }
    }
    qsort(keys, (size_t)N * N, sizeof(uint64_t), cmp_u64);

    /* First Hilbert quadrant (2 prefix bits = one 2D level). */
    prefix = keys[0];
    st = curve_key_prefix_range(keys, (size_t)N * N, 1, &prefix, 8, 2,
                                scratch_range, 2, &begin, &end);
    if (st != CURVE_OK) {
        fprintf(stderr, "prefix_range failed: %d\n", (int)st);
        return 1;
    }
    printf("Hilbert depth-1 cell of key[0] covers sorted indices [%zu, %zu)\n",
           begin, end);

    st = curve_cell_init(&cell, 2, BITS, 0);
    if (st != CURVE_OK) return 1;
    st = curve_cell_child_count(&cell, &nchild);
    printf("root child_count = %zu\n", nchild);

    for (i = 0; i < nchild; i++) {
        uint64_t root = 0;
        st = curve_cell_child(&cell, &root, (uint64_t)i, &child, child_words, 1);
        if (st != CURVE_OK) return 1;
        st = curve_hilbert_cell_bounds(&child, child_words, mins, maxs, 2,
                                       scratch, 2);
        if (st != CURVE_OK) return 1;
        printf("  child %zu: prefix=0x%016llx  AABB x[%u,%u] y[%u,%u]\n",
               i, (unsigned long long)child_words[0],
               mins[0], maxs[0], mins[1], maxs[1]);
    }
    return 0;
}
