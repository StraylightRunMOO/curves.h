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
    uint32_t hd[4], hx[4], hy[4];
    curve_morton2_batch(xs, ys, out, 4);
    curve_hilbert_xy2d_batch(xs, ys, 8, hd, 4);
    curve_hilbert_d2xy_batch(hd, 8, hx, hy, 4);

    uint32_t coords[2] = { 3u, 5u };
    uint32_t scratch[2];
    uint64_t words[1];
    const curve_status st = curve_hilbertnd_encode_words(
        coords, 2, 8, words, 1, scratch, 2);
    size_t prefix_len = 0;
    const curve_status pst = curve_key_common_prefix_bits(
        words, words, 16, &prefix_len);
    curve_cell cell;
    const curve_status cst = curve_cell_init(&cell, 2, 8, 0);
    uint64_t blo[4], bhi[4], bscratch[64];
    int bunb[4];
    uint32_t baxes[16], bmin[2] = {0u, 0u}, bmax[2] = {3u, 3u};
    size_t nint = 0;
    const curve_status rst = curve_morton_box_intervals(
        2, 4, bmin, bmax, blo, bhi, bunb, 4, &nint,
        CURVE_RANGE_EXACT, bscratch, 64, baxes, 16);

    uint32_t qbin = 99u;
    const curve_status qst = curve_quantize1_f64(
        0.5, 0.0, 1.0, 4, CURVE_QUANT_NEAREST, &qbin);
    double qv[2] = { 0.25, 0.75 }, qmin[2] = { 0.0, 0.0 }, qmax[2] = { 1.0, 1.0 };
    uint32_t qscratch[2];
    uint64_t qkey[1];
    const curve_status fst = curve_mortonnd_encode_f64(
        qv, 2, qmin, qmax, 8, CURVE_QUANT_FLOOR,
        qkey, 1, qscratch, 2, NULL, NULL, NULL);

    uint32_t px, py;
    curve_peano_d2xy(curve_peano_xy2d(1u, 2u, 1), 1, &px, &py);
    unsigned mb[2] = { 4u, 4u };
    uint64_t mkey[1];
    const curve_status mst = curve_mortonnd_encode_words_mixed(
        coords, 2, mb, mkey, 1);
    curve_dispatch disp;
    curve_dispatch_init(&disp, 0);
    uint32_t db[4];
    disp.morton2_batch(xs, ys, db, 4);

    const bool ok = (x == 3u && y == 5u
                     && out[0] == curve_morton2(1u, 4u)
                     && hd[0] == curve_hilbert_xy2d(1u, 4u, 8)
                     && CURVES_VERSION_MAJOR == 1
                     && CURVES_VERSION_MINOR == 1
                     && curve_simd_isa()[0] != '\0'
                     && curve_compiled_isa()[0] != '\0'
                     && px == 1u && py == 2u
                     && curve_gmorton2(1u, 4u) == curve_morton2(
                            1u ^ (1u >> 1), 4u ^ (4u >> 1))
                     && mst == CURVE_OK
                     && db[0] == out[0]
                     && curve_selected_isa(&disp)[0] != '\0'
                     && st == CURVE_OK
                     && pst == CURVE_OK && prefix_len == 16
                     && cst == CURVE_OK
                     && rst == CURVE_OK && nint >= 1
                     && curve_key_words_required(256, 4) == 16
                     && qst == CURVE_OK && qbin <= 15u
                     && fst == CURVE_OK);
    return ok ? 0 : 1;
}
