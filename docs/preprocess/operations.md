# Numerical operations

The preprocessing reference is implemented in ccell's `preprocess` module and
compares numerical behavior with Scanpy 1.11.5. It uses C++ native operation
contracts, not Scanpy's Python calling conventions. AnnData++ owns array
storage, views, layouts, and lifetimes. Callers supply result and scratch
buffers and own all threading. The standard reference allocates no memory,
including for sizing, validation, preparation, and empty inputs.

Every operation defines its input representation, output pointers and
capacities, scratch requirement, and write set. Checked entry points report
status, location, and required capacity without constructing exception
messages. Trusted entry points assume the documented dtype, shape, layout,
alignment, lifetime, aliasing, and sparse structure preconditions. Output
initialization and zeroing are explicit. A capacity failure does not grow or
materialize storage. Capacity and shape checks run before numerical writes when
the operation can check them at entry. A later data-dependent `QcStatus` failure
can leave output and scratch partly modified; the reference does not roll back.
Callers must use the status to determine whether outputs are valid and re-run
the operation after correcting the input if needed.

Acquisition from an AnnData++ owner may retain its existing descriptor; it must
not create numerical output or hidden workspace. Callers create arrays, attach
metadata, materialize selections, and perform any required conversion. Dense
matrix inputs may have signed strides. Sparse inputs use canonical CSR or CSC;
index and offset widths are independent, and stored zeros do not alter
mathematical values. An arbitrary gather has no direct dense pointer view and
must be materialized by the caller before an operation that needs direct
traversal.

Integer input supports summaries and count operations. Floating transforms
write FP32 or FP64 destinations. Reductions accumulate in FP64. Source and
destination overlap is operation-specific and must follow its declared
contract. The caller controls parallel decomposition, partial results, merge
order, synchronization, and workspace ownership.

## Quality control, filtering, and statistics

Quality control computes per-cell and per-feature totals and detected counts,
feature means and dropout, optional logarithmic summaries, named feature-group
totals and percentages, and optional top-feature percentages. Results remain
in input axis order. Count values use unsigned 64-bit integers; floating
summaries use FP64. QC count inputs must be finite and nonnegative.

Filtering combines minimum and maximum total-count and detected-feature bounds
into an axis-aligned mask. It does not mutate or slice the dataset. All supplied
bounds must pass.

Feature statistics return per-feature means and sample variances. An optional
observation mask limits selected rows. An empty selection produces NaNs; one
selected observation has zero variance by the reference convention.

## Normalization and transformations

Total normalization scales each nonzero row to a target using its eligible
total. The default target is the median of positive eligible row totals. An
all-zero matrix reports target one; zero-total rows remain zero. Optional
highly-expressed-feature exclusion removes a feature from size-factor
calculation when it exceeds the chosen fraction of any row's original total.
Normalization requires finite nonnegative counts. The operation reports row
factors, the included-feature mask, and the applied target into caller storage.
`normalization_factors` computes those values without applying them;
`apply_normalization` applies caller-selected row factors separately. A zero
factor leaves its row unchanged. `NormalizationOptions::trusted` skips the
per-value finite/nonnegative scan when the caller guarantees the input
contract.

Log transformation computes natural `log1p` by default, or uses the explicit
base supplied by the caller. It preserves sparse structure. Scaling computes
feature means and sample standard deviations, optionally from a row mask. It
can center, scale, and clip values. Sparse scaling without centering preserves
implicit zeros; centering requires a dense destination. Application rows and
parameter-estimation rows follow the operation's explicit mask contract.

Integer input is summarized directly. Floating transforms write FP32 or FP64;
use caller conversion when an integer source needs a floating destination.
Conversions and dense materialization are caller operations outside the
standard numerical reference.

## Highly variable features

The reference supports Seurat and Cell Ranger dispersion methods. Both compute
feature means and dispersion, group features into mean-expression bins,
normalize dispersion within each bin, and apply either cutoffs or top-score
selection. An optional batch code per row produces batch-aware selection,
selected-batch counts, and the all-batches intersection.

Seurat expects logarithmized input and reverses the specified logarithm base
before computing counts. Cell Ranger computes dispersion directly on supplied
values. Feature order determines deterministic tie results. Single-batch
top-score selection includes every feature tied at the cutoff, so the count may
exceed the requested number. Batched top-score selection ranks by selected
batch count, then normalized dispersion, then stable original feature order.
Outputs preserve original feature alignment. Intermediate and result arrays
are caller-owned.

Highly-variable-gene work can be staged: clear caller-owned batch-major moments,
accumulate row ranges into them, merge worker-local moments in caller order,
then score those moments using separately sized aligned scratch. A failed
accumulation invalidates the moments; clear and replay before merging or
scoring. The combined dense and sparse operations provide the same calculation
for one call. Score output and scratch may be partially modified on failure.

HVG scratch sizing reports a byte count. The caller supplies raw byte storage
with alignment suitable for `size_t` and `double`; scratch is transient storage,
not a live unrelated typed-object array. Input, output, and scratch storage must
not overlap for the operation's duration.

Sparse-to-sparse transforms take an explicit destination orientation. The
orientation must match the source, and the destination must carry identical
index and offset structure. Checked transforms validate the destination's
internal non-overlap and source/destination overlap; trusted calls assume these
contracts. Dense destinations may use signed strides. `log1p_into`,
`apply_normalization`, `NormalizationOptions`, and `ScaleOptions` expose a
trailing `trusted` option, defaulting to checked behavior. Use it only when the
documented storage and value preconditions are already established. The
standard operations still report capacity errors before numerical writes when
possible and do not roll back later failures.

The numerical tests and independent Scanpy fixtures cover dense, CSR, CSC,
strided inputs, FP32/FP64 calculations, integer count summaries, empty and
constant features, masks, normalization exclusions, score ordering, and ties.
See [validation](validation.md) for historical source-bound evidence and
current comparison limits.
