#pragma once

// Development-only ownership for exercising the allocation-free reference.
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace qualification {
// Byte arrays provide storage in which the kernel may create scalar objects.
struct Scratch {
    std::vector<std::byte> storage;
    void* data;
    explicit Scratch(std::size_t n) : storage(n + alignof(std::max_align_t)-1), data(nullptr) {
        auto p = reinterpret_cast<std::uintptr_t>(storage.data());
        auto aligned = (p + alignof(std::max_align_t)-1) & ~(std::uintptr_t(alignof(std::max_align_t)-1));
        data = storage.data() + (aligned-p);
    }
};
struct Qc {
    std::size_t n, g;
    std::vector<double> cell_total, gene_total, mean, dropout, log_cell_total, log_cell_detected;
    std::vector<double> log_gene_total, log_gene_mean, group_total, group_percent, group_log_total, top_percent, scratch;
    std::vector<std::uint64_t> cell_detected, gene_detected;
    std::vector<std::uint8_t> mask;
    std::size_t tops[3]{2,4,8};
    preprocess::QcOptions options;
    preprocess::QcBuffers buffers;
    Qc(anndatapp::ArrayView x, const std::vector<int>& groups) : n(x.extent(0)), g(x.extent(1)),
        cell_total(n), gene_total(g), mean(g), dropout(g), log_cell_total(n), log_cell_detected(n),
        log_gene_total(g), log_gene_mean(g), group_total(n), group_percent(n), group_log_total(n), top_percent(n*3), scratch(g),
        cell_detected(n), gene_detected(g), mask(g) {
        for (std::size_t i = 0; i < g; ++i) mask[i] = groups[i] != 0;
        options = {mask.data(), 1, tops, 3};
        buffers.cell_total=cell_total.data(); buffers.cell_detected=cell_detected.data(); buffers.cells_capacity=n;
        buffers.gene_total=gene_total.data(); buffers.gene_detected=gene_detected.data(); buffers.genes_capacity=g;
        buffers.mean=mean.data(); buffers.dropout=dropout.data(); buffers.log_cell_total=log_cell_total.data();
        buffers.log_cell_detected=log_cell_detected.data(); buffers.log_gene_total=log_gene_total.data();
        buffers.log_gene_mean=log_gene_mean.data(); buffers.group_total=group_total.data();
        buffers.group_percent=group_percent.data(); buffers.group_log_total=group_log_total.data(); buffers.group_capacity=n;
        buffers.top_percent=top_percent.data(); buffers.top_capacity=n*3;
        auto s = x.layout() == anndatapp::Layout::dense
          ? preprocess::qc(x.bytes(), x.dtype(), options, buffers, scratch.data(), scratch.size())
          : preprocess::qc(x.sparse_bytes(), x.dtype(), x.layout(), options, buffers, scratch.data(), scratch.size());
        if (!s) throw std::runtime_error("QC failed");
    }
};
struct Hvg {
    std::vector<double> mean, dispersion, normalized_dispersion;
    std::vector<std::uint8_t> selected, intersection;
    std::vector<std::size_t> nbatches, batch;
    Hvg(anndatapp::ArrayView x, const std::vector<std::int32_t>& batches, preprocess::HvgOptions options) :
      mean(x.extent(1)), dispersion(x.extent(1)), normalized_dispersion(x.extent(1)),
      selected(x.extent(1)), intersection(x.extent(1)), nbatches(x.extent(1)), batch(x.extent(0)) {
        std::vector<std::int32_t> codes;
        if (options.batched) for (std::size_t i=0; i<batch.size(); ++i) {
            std::size_t j=0; for (; j<codes.size() && codes[j]!=batches[i]; ++j) {}
            if (j==codes.size()) codes.push_back(batches[i]); batch[i]=j;
        }
        std::size_t count=options.batched ? codes.size() : 1, bytes=0;
        auto s=preprocess::hvg_scratch_bytes(x.extent(1),count,options,bytes);
        if (s.code!=preprocess::HvgCode::ok) throw std::runtime_error("HVG sizing failed");
        Scratch scratch(bytes);
        preprocess::HvgOutput out{mean.data(),dispersion.data(),normalized_dispersion.data(),selected.data(),intersection.data(),nbatches.data(),x.extent(1)};
        const auto* b=options.batched ? batch.data() : nullptr;
        if (x.layout()==anndatapp::Layout::dense) {
            auto a=x.bytes(); s=preprocess::hvg_dense(a.data(),x.dtype(),x.extent(0),x.extent(1),a.stride_bytes(0),a.stride_bytes(1),b,count,options,out,scratch.data,bytes);
        } else {
            auto a=x.sparse_bytes(); s=preprocess::hvg_sparse(a.values(),x.dtype(),a.indices(),a.offsets(),a.index_size(),a.offset_size(),x.layout(),x.extent(0),x.extent(1),a.nnz(),b,count,options,out,scratch.data,bytes);
        }
        if (s.code!=preprocess::HvgCode::ok) throw std::runtime_error("HVG calculation failed");
    }
};
} // namespace qualification
