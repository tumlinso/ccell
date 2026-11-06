module;
#include <stddef.h>
#include <stdint.h>
export module preprocess;
export import preprocess.qc;
export import preprocess.hvg;

export namespace preprocess {

// Compose arbitrary caller state without imposing reference layouts. A zero
// block size reports invalid_argument; bindings and synchronization are caller-owned.
template<class Operation> QcStatus for_blocks(size_t total, size_t block, Operation& operation) {
    if (!block) return {QcCode::invalid_argument};
    for (size_t first = 0; first < total;) {
        size_t count = total-first < block ? total-first : block;
        operation(first,count); first += count;
    }
    return {};
}
QcStatus for_blocks(size_t total, size_t block, void (*operation)(void*, size_t, size_t), void* context);

// Caller acquires valid rank-2 AnnData byte descriptors before entering these
// routines. Dense byte strides must be representable. Sparse inputs must have
// canonical sorted-unique indices and valid offsets. Buffers must not overlap,
// except identical source/destination storage with matching shape and strides.
// Acquired descriptor lifetime and concurrent access are caller responsibilities.
struct NormalizationOptions {
    bool target_is_set = false;
    double target_sum = 1;
    bool exclude_highly_expressed = false;
    double max_fraction = .05;
    bool trusted = false;
};
// Required capacities are rows factors, columns included flags, rows scratch.
// Inputs and buffers are borrowed; factor and included contents are overwritten.
// trusted=true omits per-element finite/nonnegative checks. Caller guarantees
// finite nonnegative values, representable finite row sums and valid structure.
struct NormalizationBuffers {
    double* factors = nullptr; size_t factors_capacity = 0;
    uint8_t* included = nullptr; size_t included_capacity = 0;
    double* scratch = nullptr; size_t scratch_capacity = 0;
    double target_sum = 1;
};
using DenseOutput = decltype(((const anndatapp::ArrayEdit*)nullptr)->bytes());
using SparseOutput = decltype(((const anndatapp::ArrayEdit*)nullptr)->sparse_bytes());
// Sparse output must declare the same orientation and contain exactly the input
// index/offset structure. Checked calls compare structures; trusted calls assume it.

QcStatus normalization_factors(const DenseInput&, anndatapp::DType, const NormalizationOptions&, NormalizationBuffers&);
QcStatus normalization_factors(const SparseInput&, anndatapp::DType, anndatapp::Layout, const NormalizationOptions&, NormalizationBuffers&);
QcStatus normalize_into(const DenseInput&, anndatapp::DType, const DenseOutput&, anndatapp::DType, const NormalizationOptions&, NormalizationBuffers&);
QcStatus normalize_into(const SparseInput&, anndatapp::DType, anndatapp::Layout, const DenseOutput&, anndatapp::DType, const NormalizationOptions&, NormalizationBuffers&);
QcStatus normalize_into(const SparseInput&, anndatapp::DType, anndatapp::Layout, const SparseOutput&, anndatapp::DType, anndatapp::Layout destination_layout, const NormalizationOptions&, NormalizationBuffers&);

// Apply caller-supplied per-row factors; a zero factor leaves that row unchanged.
// trusted=true skips the finite/nonnegative factor scan. Input and output contracts
// are the same as normalize_into. The factors array has at least rows elements.
QcStatus apply_normalization(const DenseInput&, anndatapp::DType, const DenseOutput&, anndatapp::DType, const double* factors, size_t capacity, bool trusted = false);
QcStatus apply_normalization(const SparseInput&, anndatapp::DType, anndatapp::Layout, const DenseOutput&, anndatapp::DType, const double* factors, size_t capacity, bool trusted = false);
QcStatus apply_normalization(const SparseInput&, anndatapp::DType, anndatapp::Layout, const SparseOutput&, anndatapp::DType, anndatapp::Layout destination_layout, const double* factors, size_t capacity, bool trusted = false);

// base=0 specifies the natural logarithm; other bases must be positive and != 1.
// trusted=true assumes the destination has valid nonoverlapping elements.
QcStatus log1p_into(const DenseInput&, anndatapp::DType, const DenseOutput&, anndatapp::DType, double base = 0, bool trusted = false);
QcStatus log1p_into(const SparseInput&, anndatapp::DType, anndatapp::Layout, const DenseOutput&, anndatapp::DType, double base = 0, bool trusted = false);
QcStatus log1p_into(const SparseInput&, anndatapp::DType, anndatapp::Layout, const SparseOutput&, anndatapp::DType, anndatapp::Layout destination_layout, double base = 0, bool trusted = false);

struct ScaleOptions {
    bool zero_center = true, clip = false;
    double max_value = 0;
    const uint8_t* observations = nullptr;
    bool trusted = false;
};
// Fit requires columns elements in mean/stddev and uint64 scratch. Input mask
// selects rows for both estimation and application. Sparse centered output is dense.
// trusted=true omits per-element finite/structure validation in statistics.
// Caller guarantees finite input and canonical structure, valid mask and parameters.
QcStatus scale_into(const DenseInput&, anndatapp::DType, const DenseOutput&, anndatapp::DType, const ScaleOptions&, double* mean, double* stddev, size_t capacity, uint64_t* scratch, size_t scratch_capacity);
QcStatus scale_into(const SparseInput&, anndatapp::DType, anndatapp::Layout, const DenseOutput&, anndatapp::DType, const ScaleOptions&, double* mean, double* stddev, size_t capacity, uint64_t* scratch, size_t scratch_capacity);
QcStatus scale_into(const SparseInput&, anndatapp::DType, anndatapp::Layout, const SparseOutput&, anndatapp::DType, anndatapp::Layout destination_layout, const ScaleOptions&, double* mean, double* stddev, size_t capacity, uint64_t* scratch, size_t scratch_capacity);
QcStatus apply_scale(const DenseInput&, anndatapp::DType, const DenseOutput&, anndatapp::DType, const ScaleOptions&, const double* mean, const double* stddev, size_t capacity);
QcStatus apply_scale(const SparseInput&, anndatapp::DType, anndatapp::Layout, const DenseOutput&, anndatapp::DType, const ScaleOptions&, const double* mean, const double* stddev, size_t capacity);
QcStatus apply_scale(const SparseInput&, anndatapp::DType, anndatapp::Layout, const SparseOutput&, anndatapp::DType, anndatapp::Layout destination_layout, const ScaleOptions&, const double* mean, const double* stddev, size_t capacity);
}
