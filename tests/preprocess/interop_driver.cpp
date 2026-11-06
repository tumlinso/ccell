#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

import preprocess;

#include "qualification_helpers.hpp"

namespace {

template<class T> void emit_values(anndatapp::ArrayView a) {
    auto v = a.dense<T>();
    if (v.rank() == 1) {
        for (std::size_t i = 0; i < v.extent(0); ++i) std::cout << ' ' << double(v.data()[i*v.stride(0)]);
    } else if (v.rank() == 2) {
        for (std::size_t i = 0; i < v.extent(0); ++i)
            for (std::size_t j = 0; j < v.extent(1); ++j)
                std::cout << ' ' << double(v.data()[i*v.stride(0) + j*v.stride(1)]);
    } else throw std::runtime_error("unsupported output rank");
}

void emit(const char* name, anndatapp::ArrayView a) {
    std::cout << name << ' ' << a.rank();
    for (std::size_t i = 0; i < a.rank(); ++i) std::cout << ' ' << a.extent(i);
    switch (a.dtype()) {
    case anndatapp::DType::boolean: emit_values<bool>(a); break;
    case anndatapp::DType::int8: emit_values<std::int8_t>(a); break;
    case anndatapp::DType::uint8: emit_values<std::uint8_t>(a); break;
    case anndatapp::DType::int32: emit_values<std::int32_t>(a); break;
    case anndatapp::DType::uint32: emit_values<std::uint32_t>(a); break;
    case anndatapp::DType::int64: emit_values<std::int64_t>(a); break;
    case anndatapp::DType::uint64: emit_values<std::uint64_t>(a); break;
    case anndatapp::DType::float32: emit_values<float>(a); break;
    case anndatapp::DType::float64: emit_values<double>(a); break;
    default: throw std::runtime_error("unsupported output dtype");
    }
    std::cout << '\n';
}

void emit(const char* name, const anndatapp::Array& a) { emit(name, a.view()); }

void emit_matrix(anndatapp::ArrayView a) {
    const auto n = a.extent(0), g = a.extent(1);
    std::vector<double> values(n*g, 0);
    auto read = [&]<class T>() {
        if (a.layout() == anndatapp::Layout::dense) {
            auto v = a.dense<T>();
            for (std::size_t i = 0; i < n; ++i) for (std::size_t j = 0; j < g; ++j)
                values[i*g+j] = double(v.data()[std::ptrdiff_t(i)*v.stride(0)+std::ptrdiff_t(j)*v.stride(1)]);
        } else if (a.layout() == anndatapp::Layout::csr) {
            auto v = a.csr<T, std::int32_t, std::int64_t>();
            for (std::size_t i = 0; i < n; ++i)
                for (auto k = v.offsets()[i]; k < v.offsets()[i+1]; ++k)
                    values[i*g+v.indices()[k]] = double(v.values()[k]);
        } else {
            auto v = a.csc<T, std::int32_t, std::int64_t>();
            for (std::size_t j = 0; j < g; ++j)
                for (auto k = v.offsets()[j]; k < v.offsets()[j+1]; ++k)
                    values[v.indices()[k]*g+j] = double(v.values()[k]);
        }
    };
    switch (a.dtype()) {
    case anndatapp::DType::int32: read.template operator()<std::int32_t>(); break;
    case anndatapp::DType::float32: read.template operator()<float>(); break;
    case anndatapp::DType::float64: read.template operator()<double>(); break;
    default: throw std::runtime_error("unsupported matrix dtype");
    }
    std::cout << "matrix 2 " << n << ' ' << g;
    for (double v : values) std::cout << ' ' << v;
    std::cout << '\n';
}

template<class T> void emit(const char* name, const std::vector<T>& values) {
    std::cout << name << " 1 " << values.size();
    for (auto v : values) std::cout << ' ' << double(v);
    std::cout << '\n';
}
void emit_top(const char* name, const std::vector<double>& values, std::size_t rows, std::size_t tops) {
    std::cout << name << " 2 " << rows << ' ' << tops;
    for (auto v : values) std::cout << ' ' << v;
    std::cout << '\n';
}

#include "annotation_fixture.hpp"

template<class T>
void run(const std::string& operation, int layout, int mode, std::size_t n, std::size_t g,
         const std::vector<double>& input, const std::vector<int>& group_values,
         const std::vector<int>& observation_values, const std::vector<std::int32_t>& batch_values) {
    std::vector<T> values;
    std::vector<std::int32_t> indices;
    std::vector<std::int64_t> offsets;
    anndatapp::Array a;
    if (layout == 0 || layout == 3) {
        values.resize(n*g);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < g; ++j)
                values[(layout == 3 ? n - 1 - i : i)*g + j] = T(input[i*g + j]);
        T* p = values.data();
        if (layout == 3 && n) p += (n-1)*g;
        a = anndatapp::Array::borrow_dense(p, {n,g}, {layout == 3 ? -std::ptrdiff_t(g) : std::ptrdiff_t(g), 1});
    } else {
        const bool csr = layout == 1;
        const std::size_t major = csr ? n : g, minor = csr ? g : n;
        offsets.push_back(0);
        for (std::size_t i = 0; i < major; ++i) {
            for (std::size_t j = 0; j < minor; ++j) {
                const std::size_t row = csr ? i : j, col = csr ? j : i;
                T x = T(input[row*g + col]);
                if (x != T(0) || (row*g + col) % 31 == 0) {
                    values.push_back(x); indices.push_back(std::int32_t(j));
                }
            }
            offsets.push_back(std::int64_t(values.size()));
        }
        a = anndatapp::Array::borrow_sparse(csr ? anndatapp::Layout::csr : anndatapp::Layout::csc,
            n, g, values.size(), values.data(), indices.data(), offsets.data());
    }
    if (operation == "annotations") {
        if (layout != 0 || std::is_integral_v<T>) throw std::invalid_argument("annotations require dense floating input");
        run_annotations(a, mode, n, g, group_values);
        return;
    }
    auto x = a.view();
    std::vector<std::uint8_t> mask(n);
    for (std::size_t i=0; i<n; ++i) mask[i] = observation_values[i] != 0;
    if (operation == "qc" || operation == "filters") {
        qualification::Qc q(x, group_values);
        if (operation == "filters") {
            preprocess::FilterBounds cells{true,true,true,true,4,10000,2,g};
            preprocess::FilterBounds genes{true,true,true,true,5,10000,2,n};
            std::vector<std::uint8_t> c(n), v(g);
            if (!preprocess::filter(q.cell_total.data(),q.cell_detected.data(),n,cells,c.data(),c.size()) ||
                !preprocess::filter(q.gene_total.data(),q.gene_detected.data(),g,genes,v.data(),v.size()))
                throw std::runtime_error("filter failed");
            emit("cells",c); emit("genes",v); return;
        }
        emit("cell_total",q.cell_total); emit("cell_detected",q.cell_detected);
        emit("gene_total",q.gene_total); emit("gene_detected",q.gene_detected);
        emit("gene_mean",q.mean); emit("gene_dropout",q.dropout);
        emit("log_cell_total",q.log_cell_total); emit("log_cell_detected",q.log_cell_detected);
        emit("log_gene_total",q.log_gene_total); emit("log_gene_mean",q.log_gene_mean);
        emit("group_total",q.group_total); emit("group_percent",q.group_percent); emit("group_log_total",q.group_log_total);
        emit_top("top_percent",q.top_percent,n,3);
    } else if (operation == "stats") {
        std::vector<double> mean(g), variance(g); std::vector<std::uint64_t> scratch(g); std::size_t observations=0;
        const auto* m = mode ? mask.data() : nullptr;
        auto status = x.layout() == anndatapp::Layout::dense
          ? preprocess::feature_stats(x.bytes(),x.dtype(),m,mean.data(),variance.data(),g,scratch.data(),scratch.size(),&observations)
          : preprocess::feature_stats(x.sparse_bytes(),x.dtype(),x.layout(),m,mean.data(),variance.data(),g,scratch.data(),scratch.size(),&observations);
        if (!status) throw std::runtime_error("feature_stats failed");
        emit("mean",mean); emit("variance",variance);
    } else if (operation == "hvg") {
        preprocess::HvgOptions o;
        o.method=(mode & 1) ? preprocess::HvgMethod::cell_ranger : preprocess::HvgMethod::seurat;
        o.n_bins=10; o.has_top=!(mode & 4); o.n_top_genes=12; o.batched=mode & 2;
        o.has_log_base=mode & (8|16); o.log_base=(mode & 8) ? 2 : .5;
        qualification::Hvg h(x,batch_values,o);
        emit("mean",h.mean); emit("dispersion",h.dispersion); emit("normalized_dispersion",h.normalized_dispersion);
        emit("selected",h.selected); if (o.batched) { emit("nbatches",h.nbatches); emit("intersection",h.intersection); }
    } else if (operation == "normalize" || operation == "normalize_into" || operation == "log1p" || operation == "log1p_into" || operation == "scale") {
        const bool single = operation == "normalize_into" || operation == "log1p_into" || (std::is_same_v<T,float> && !(operation == "scale" && layout != 0 && layout != 3));
        anndatapp::Array output = single ? anndatapp::Array::empty<float>({n,g}) : anndatapp::Array::empty<double>({n,g});
        auto y=output.edit(); auto dst=y.bytes(); const auto dt=output.dtype(); preprocess::QcStatus status;
        if (operation == "normalize" || operation == "normalize_into") {
            preprocess::NormalizationOptions o; o.target_is_set=mode!=1; o.target_sum=100;
            o.exclude_highly_expressed=mode==2; o.max_fraction=.2;
            std::vector<double> factors(n), scratch(n); std::vector<std::uint8_t> included(g);
            preprocess::NormalizationBuffers b{factors.data(),n,included.data(),g,scratch.data(),n};
            status=x.layout()==anndatapp::Layout::dense ? preprocess::normalize_into(x.bytes(),x.dtype(),dst,dt,o,b)
              : preprocess::normalize_into(x.sparse_bytes(),x.dtype(),x.layout(),dst,dt,o,b);
            if (!status) throw std::runtime_error("normalize failed");
            emit_matrix(output.view()); emit("factors",factors); emit("included",included);
            std::cout << "target 1 1 " << b.target_sum << '\n';
        } else if (operation == "log1p" || operation == "log1p_into") {
            double base=mode ? (mode==1 ? 2.0 : .5) : 0;
            status=x.layout()==anndatapp::Layout::dense ? preprocess::log1p_into(x.bytes(),x.dtype(),dst,dt,base)
              : preprocess::log1p_into(x.sparse_bytes(),x.dtype(),x.layout(),dst,dt,base);
            if (!status) throw std::runtime_error("log1p failed"); emit_matrix(output.view());
        } else {
            preprocess::ScaleOptions o; o.zero_center=mode!=1; o.clip=mode!=0; o.max_value=3;
            o.observations=mode==2 ? mask.data() : nullptr;
            std::vector<double> mean(g), dev(g); std::vector<std::uint64_t> scratch(g);
            status=x.layout()==anndatapp::Layout::dense ? preprocess::scale_into(x.bytes(),x.dtype(),dst,dt,o,mean.data(),dev.data(),g,scratch.data(),g)
              : preprocess::scale_into(x.sparse_bytes(),x.dtype(),x.layout(),dst,dt,o,mean.data(),dev.data(),g,scratch.data(),g);
            if (!status) throw std::runtime_error("scale failed"); emit_matrix(output.view()); emit("mean",mean); emit("stddev",dev);
        }
    } else throw std::invalid_argument("unknown operation");

}

} // namespace

int main(int argc, char** argv) try {
    if (argc != 2) throw std::invalid_argument("usage: ccell_preprocess_interop_driver operation < input");
    std::size_t n, g; int dtype, layout, mode;
    if (!(std::cin >> n >> g >> dtype >> layout >> mode) || n > 10000 || g > 10000)
        throw std::invalid_argument("invalid header");
    std::vector<double> values(n*g);
    std::vector<int> groups(g), observations(n);
    std::vector<std::int32_t> batches(n);
    for (auto& value : values) std::cin >> value;
    for (auto& value : groups) std::cin >> value;
    for (auto& value : observations) std::cin >> value;
    for (auto& value : batches) std::cin >> value;
    if (!std::cin) throw std::invalid_argument("truncated input");
    std::cout << std::setprecision(17);
    if (dtype == 32) run<float>(argv[1], layout, mode, n, g, values, groups, observations, batches);
    else if (dtype == 64) run<double>(argv[1], layout, mode, n, g, values, groups, observations, batches);
    else if (dtype == 0) run<std::int32_t>(argv[1], layout, mode, n, g, values, groups, observations, batches);
    else throw std::invalid_argument("unsupported dtype");
    return 0;
} catch (const std::exception& e) {
    std::cerr << "interop: " << e.what() << '\n';
    return 1;
}
