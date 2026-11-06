#include <stddef.h>
#include <stdint.h>

import anndatapp.core;
import preprocess;

int main() {
    auto a = anndatapp::Array::empty<float>({3,3});
    auto x = a.edit().dense<float>(); auto* p = x.data();
    p[0] = 2; p[1] = 0; p[2] = 1;
    p[3] = 4; p[4] = 2; p[5] = 0;
    p[6] = 0; p[7] = 1; p[8] = 3;

    auto v = a.view(); auto d = v.bytes();
    double cell_total[3], gene_total[3], mean[3], dropout[3];
    double log_cell_total[3], log_cell_detected[3], log_gene_total[3], log_gene_mean[3];
    uint64_t cell_detected[3], gene_detected[3];
    double scratch[3];
    preprocess::QcBuffers out{
        .cell_total = cell_total, .cell_detected = cell_detected, .cells_capacity = 3,
        .gene_total = gene_total, .gene_detected = gene_detected, .genes_capacity = 3,
        .mean = mean, .dropout = dropout, .log_cell_total = log_cell_total,
        .log_cell_detected = log_cell_detected, .log_gene_total = log_gene_total,
        .log_gene_mean = log_gene_mean};
    auto status = preprocess::qc(d, anndatapp::DType::float32, {}, out, scratch, 3);
    if (!status) return 1;

    double variance[3]; uint64_t stored[3]; size_t observations;
    status = preprocess::feature_stats(d, anndatapp::DType::float32, nullptr,
                                       mean, variance, 3, stored, 3, &observations);
    if (!status) return 2;
    if (observations != 3 || cell_total[0] != 3 || cell_total[1] != 6 ||
        gene_total[0] != 6 || mean[0] != 2 || variance[0] != 4) return 3;

    uint8_t keep[3];
    status = preprocess::filter(cell_total, cell_detected, 3,
                                {.has_min_counts = true, .min_counts = 4}, keep, 3);
    if (!status) return 4;
    return keep[0] == 0 && keep[1] == 1 && keep[2] == 1 ? 0 : 5;
}
