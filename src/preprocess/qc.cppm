module;
#include <stddef.h>
#include <stdint.h>
export module preprocess.qc;
export import anndatapp.core;
export namespace preprocess {

using DenseInput = decltype(((const anndatapp::ArrayView*)nullptr)->bytes());
using SparseInput = decltype(((const anndatapp::ArrayView*)nullptr)->sparse_bytes());

// The caller acquires descriptors before entering reference code. Dense inputs
// must have rank two and byte strides representable as ptrdiff_t; malformed
// AnnData descriptor acquisition/stride failures are caller responsibilities.

enum class QcCode : uint8_t { ok, invalid_argument, insufficient_capacity, invalid_value, invalid_structure, unsupported_type };
struct QcStatus { QcCode code = QcCode::ok; size_t index = 0, required = 0; explicit operator bool() const { return code == QcCode::ok; } };

// Group masks are contiguous group-major boolean bytes with `columns` bytes per group.
struct QcOptions { const uint8_t* masks = nullptr; size_t groups = 0; const size_t* percent_top = nullptr; size_t tops = 0; };
// All arrays are caller-owned; optional arrays may be null. Cell arrays have rows
// entries, gene arrays columns entries; group arrays use group-major rows;
// top_percent uses row-major tops. Capacities are element counts.
struct QcBuffers {
    double* cell_total = nullptr; uint64_t* cell_detected = nullptr; size_t cells_capacity = 0;
    double* gene_total = nullptr; uint64_t* gene_detected = nullptr; size_t genes_capacity = 0;
    double* mean = nullptr; double* dropout = nullptr;
    double* log_cell_total = nullptr; double* log_cell_detected = nullptr;
    double* log_gene_total = nullptr; double* log_gene_mean = nullptr;
    double* group_total = nullptr; double* group_percent = nullptr; double* group_log_total = nullptr; size_t group_capacity = 0;
    double* top_percent = nullptr; size_t top_capacity = 0;
};
// Input descriptors are acquired by caller. Partial range [first,last) accumulates
// into caller-zeroed buffers and writes only its cell rows; invalid input may leave
// partial values. Finalize and merge are separate and perform no allocation.
QcStatus qc_check(size_t rows, size_t columns, const QcOptions&, const QcBuffers&, size_t scratch_elements);
// Validate finite values and canonical sparse structure without allocation.
QcStatus qc_validate(const DenseInput&, anndatapp::DType, bool nonnegative);
QcStatus qc_validate(const SparseInput&, anndatapp::DType, anndatapp::Layout, bool nonnegative);
void qc_zero(size_t rows, size_t columns, size_t groups, const QcBuffers&);
QcStatus qc_accumulate(const DenseInput&, anndatapp::DType, size_t first, size_t last, const QcOptions&, const QcBuffers&, double* scratch, size_t scratch_elements);
QcStatus qc_accumulate(const SparseInput&, anndatapp::DType, anndatapp::Layout, size_t first, size_t last, const QcOptions&, const QcBuffers&, double* scratch, size_t scratch_elements);
// Caller guarantees valid finite nonnegative data, canonical sparse structure,
// valid outputs and ranges; no per-entry validation occurs.
QcStatus qc_accumulate_unchecked(const DenseInput&, anndatapp::DType, size_t first, size_t last, const QcOptions&, const QcBuffers&, double* scratch);
QcStatus qc_accumulate_unchecked(const SparseInput&, anndatapp::DType, anndatapp::Layout, size_t first, size_t last, const QcOptions&, const QcBuffers&, double* scratch);
void qc_merge(size_t columns, double* total, uint64_t* detected, const double* part_total, const uint64_t* part_detected);
void qc_finish(size_t rows, size_t columns, const QcOptions&, const QcBuffers&);
QcStatus qc(const DenseInput&, anndatapp::DType, const QcOptions&, const QcBuffers&, double* scratch, size_t scratch_elements);
QcStatus qc(const SparseInput&, anndatapp::DType, anndatapp::Layout, const QcOptions&, const QcBuffers&, double* scratch, size_t scratch_elements);

struct FilterBounds { bool has_min_counts = false, has_max_counts = false, has_min_detected = false, has_max_detected = false; double min_counts = 0, max_counts = 0; uint64_t min_detected = 0, max_detected = 0; };
QcStatus filter(const double* totals, const uint64_t* detected, size_t count, const FilterBounds&, uint8_t* keep, size_t capacity);
// Welford partials represent selected rows in [first,last). Each partial owns
// no memory: caller supplies per-gene mean/M2 arrays and count scratch. Separate
// workers use disjoint arrays; merge in the caller's chosen order.
struct StatsPartial { size_t observations = 0, columns = 0; double* mean = nullptr; double* m2 = nullptr; size_t capacity = 0; };
QcStatus stats_partial(const DenseInput&, anndatapp::DType, size_t first, size_t last, const uint8_t* mask, StatsPartial&, uint64_t* scratch, size_t scratch_elements);
QcStatus stats_partial(const SparseInput&, anndatapp::DType, anndatapp::Layout, size_t first, size_t last, const uint8_t* mask, StatsPartial&, uint64_t* scratch, size_t scratch_elements);
QcStatus stats_merge(StatsPartial&, const StatsPartial&);
// Finalizing an empty partial sets its mean and the output variance to NaN.
QcStatus stats_finish(StatsPartial&, double* variance, size_t capacity);
// Scratch holds `columns` uint64_t values; mask is optional contiguous rows bytes.
// Inputs must be canonical sparse with sorted unique indices. All values must be finite.
QcStatus feature_stats(const DenseInput&, anndatapp::DType, const uint8_t* mask, double* mean, double* variance, size_t capacity, uint64_t* scratch, size_t scratch_elements, size_t* observations);
QcStatus feature_stats(const SparseInput&, anndatapp::DType, anndatapp::Layout, const uint8_t* mask, double* mean, double* variance, size_t capacity, uint64_t* scratch, size_t scratch_elements, size_t* observations);
// Assumes finite values and canonical sparse structure; capacity checks remain.
QcStatus feature_stats_unchecked(const DenseInput&, anndatapp::DType, const uint8_t* mask, double* mean, double* variance, size_t capacity, uint64_t* scratch, size_t scratch_elements, size_t* observations);
QcStatus feature_stats_unchecked(const SparseInput&, anndatapp::DType, anndatapp::Layout, const uint8_t* mask, double* mean, double* variance, size_t capacity, uint64_t* scratch, size_t scratch_elements, size_t* observations);

}
