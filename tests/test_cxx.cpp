/* Compile-only + smoke test that curves.h is C++11-safe. */
#include "curves.h"

int main()
{
    uint32_t x = 0, y = 0;
    uint32_t d = curve_hilbert_xy2d(3, 5, 8);
    curve_hilbert_d2xy(d, 8, &x, &y);

    uint32_t xs[4] = {1u, 2u, 3u, 4u};
    uint32_t ys[4] = {4u, 3u, 2u, 1u};
    uint32_t out[4];
    curve_morton2_batch(xs, ys, out, 4);

    const bool ok = (x == 3u && y == 5u
                     && out[0] == curve_morton2(1u, 4u)
                     && CURVES_VERSION_MAJOR == 1
                     && curve_simd_isa()[0] != '\0');
    return ok ? 0 : 1;
}
