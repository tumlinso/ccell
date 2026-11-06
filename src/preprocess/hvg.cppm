module;
#include <stddef.h>
#include <stdint.h>

export module preprocess.hvg;
export import anndatapp.core;

export namespace preprocess {

enum class HvgMethod { seurat, cell_ranger };
enum class HvgCode { ok, invalid_argument, insufficient_capacity, overflow, invalid_input, empty_batch, no_scores, duplicate_bin_edges };
struct HvgStatus { HvgCode code = HvgCode::ok; size_t index = 0, required_bytes = 0; };
struct HvgOptions {
    HvgMethod method = HvgMethod::seurat;
    size_t n_top_genes = 0, n_bins = 20;
    double min_mean = .0125, max_mean = 3, min_disp = .5, max_disp = __builtin_inf();
    double log_base = 0;
    bool has_top = false, has_log_base = false, batched = false;
};
struct HvgOutput {
    double *mean = nullptr, *dispersion = nullptr, *normalized_dispersion = nullptr;
    uint8_t *selected = nullptr, *intersection = nullptr;
    size_t* nbatches = nullptr;
    size_t capacity = 0;
};
// Batch-major moment arrays have capacity >= genes*batches; counts has
// capacity >= batches. Arrays are owned by the caller and must not overlap.
struct HvgMoments {
    double *sum = nullptr, *square = nullptr;
    uint8_t* expressed = nullptr;
    size_t* counts = nullptr;
    size_t capacity = 0, batch_capacity = 0;
};

// Scratch is suitably aligned raw or byte storage permitting implicit-lifetime
// double and size_t objects, aligned to alignof(size_t) and alignof(double).
// Batch codes are contiguous [0,batches), caller-chosen; every batch is used.
// Input, output, and scratch storage must not overlap; all pointers remain valid
// until return. Capacities are in elements for output, bytes for scratch.
// On insufficient capacity no output is written; other failures may leave
// partial output and scratch contents. Batched calls also fill nbatches and
// intersection; unbatched calls leave those two arrays untouched.
HvgStatus hvg_scratch_bytes(size_t genes, size_t batches, const HvgOptions&, size_t& bytes);
// Rank batched genes by descending selected-batch count, then descending
// normalized dispersion, then original gene position. order has capacity genes.
HvgStatus hvg_rank_batched(const double* normalized_dispersion, const size_t* nbatches,
    size_t genes, size_t n_top_genes, uint8_t* selected, size_t* order, size_t capacity);
// Clear before the first accumulation; successive calls add rows. Input row
// batch codes address this same batch space. Any failed accumulation invalidates
// all moments: clear and replay before further accumulation, merge, or score.
HvgStatus hvg_moments_clear(size_t genes, size_t batches, HvgMoments);
// Add source into destination. Both refer to the same gene/batch order; arrays
// may be identical pairwise or disjoint, and counts must not overflow.
HvgStatus hvg_moments_merge(size_t genes, size_t batches, HvgMoments destination, const HvgMoments& source);
HvgStatus hvg_accumulate_dense(const void* values, anndatapp::DType, size_t rows, size_t genes,
    ptrdiff_t row_stride_bytes, ptrdiff_t gene_stride_bytes, const size_t* batch, size_t batches,
    const HvgOptions&, HvgMoments);
HvgStatus hvg_accumulate_dense_unchecked(const void* values, anndatapp::DType, size_t rows, size_t genes,
    ptrdiff_t row_stride_bytes, ptrdiff_t gene_stride_bytes, const size_t* batch, size_t batches,
    const HvgOptions&, HvgMoments);
HvgStatus hvg_accumulate_sparse(const void* values, anndatapp::DType, const void* indices, const void* offsets,
    size_t index_width, size_t offset_width, anndatapp::Layout, size_t rows, size_t genes,
    size_t nnz, const size_t* batch, size_t batches, const HvgOptions&, HvgMoments);
HvgStatus hvg_accumulate_sparse_unchecked(const void* values, anndatapp::DType, const void* indices, const void* offsets,
    size_t index_width, size_t offset_width, anndatapp::Layout, size_t rows, size_t genes,
    size_t nnz, const size_t* batch, size_t batches, const HvgOptions&, HvgMoments);
// Score caller-constructed moments without copying them. Input/scratch/output
// must be disjoint. Each batch must have nonzero count. Use the scoring scratch
// query; it excludes input moments and requires the alignment documented above.
HvgStatus hvg_score_scratch_bytes(size_t genes, const HvgOptions&, size_t& bytes);
HvgStatus hvg_score(const HvgMoments&, size_t genes, size_t batches, const HvgOptions&,
    HvgOutput, void* scratch, size_t scratch_bytes);
HvgStatus hvg_score_unchecked(const HvgMoments&, size_t genes, size_t batches, const HvgOptions&,
    HvgOutput, void* scratch, size_t scratch_bytes);
HvgStatus hvg_dense(const void* values, anndatapp::DType, size_t rows, size_t genes,
    ptrdiff_t row_stride_bytes, ptrdiff_t gene_stride_bytes, const size_t* batch, size_t batches,
    const HvgOptions&, HvgOutput, void* scratch, size_t scratch_bytes);
HvgStatus hvg_sparse(const void* values, anndatapp::DType, const void* indices, const void* offsets,
    size_t index_width, size_t offset_width, anndatapp::Layout, size_t rows, size_t genes,
    size_t nnz, const size_t* batch, size_t batches, const HvgOptions&, HvgOutput,
    void* scratch, size_t scratch_bytes);

// Trusted variants omit value and sparse-structure checks. The caller supplies
// finite, nonnegative transformed values and canonical compressed structures.
HvgStatus hvg_dense_unchecked(const void* values, anndatapp::DType, size_t rows, size_t genes,
    ptrdiff_t row_stride_bytes, ptrdiff_t gene_stride_bytes, const size_t* batch, size_t batches,
    const HvgOptions&, HvgOutput, void* scratch, size_t scratch_bytes);
HvgStatus hvg_sparse_unchecked(const void* values, anndatapp::DType, const void* indices, const void* offsets,
    size_t index_width, size_t offset_width, anndatapp::Layout, size_t rows, size_t genes,
    size_t nnz, const size_t* batch, size_t batches, const HvgOptions&, HvgOutput,
    void* scratch, size_t scratch_bytes);

} // namespace preprocess
