module;
#include <stddef.h>
#include <stdint.h>
module preprocess;

namespace preprocess {

QcStatus for_blocks(size_t total, size_t block, void (*operation)(void*, size_t, size_t), void* context) {
    if (!block || !operation) return {QcCode::invalid_argument};
    for (size_t first = 0; first < total;) {
        size_t count = total-first < block ? total-first : block;
        operation(context,first,count); first += count;
    }
    return {};
}

namespace {
template<class T, class Byte> T load(const Byte* p) { T x; __builtin_memcpy(&x,p,sizeof(T)); return x; }
template<class T, class Byte> void store(Byte* p, T x) { __builtin_memcpy(p,&x,sizeof(T)); }
size_t width(anndatapp::DType t) {
    using D = anndatapp::DType;
    switch (t) {
    case D::boolean: case D::int8: case D::uint8: return 1;
    case D::int16: case D::uint16: return 2;
    case D::int32: case D::uint32: case D::float32: return 4;
    case D::int64: case D::uint64: case D::float64: return 8;
    }
    return 0;
}
template<class F> QcStatus dispatch(anndatapp::DType type, F f) {
    using D = anndatapp::DType;
    switch (type) {
    case D::boolean: return f.template operator()<bool>();
    case D::int8: return f.template operator()<int8_t>();
    case D::uint8: return f.template operator()<uint8_t>();
    case D::int16: return f.template operator()<int16_t>();
    case D::uint16: return f.template operator()<uint16_t>();
    case D::int32: return f.template operator()<int32_t>();
    case D::uint32: return f.template operator()<uint32_t>();
    case D::int64: return f.template operator()<int64_t>();
    case D::uint64: return f.template operator()<uint64_t>();
    case D::float32: return f.template operator()<float>();
    case D::float64: return f.template operator()<double>();
    }
    return {QcCode::unsupported_type};
}
template<class T, class F> void entries(const DenseInput& x, F f) {
    size_t m = x.extent(0), n = x.extent(1); auto* p = x.data();
    ptrdiff_t rs = x.stride_bytes(0), cs = x.stride_bytes(1);
    for (size_t i = 0; i < m; ++i) for (size_t j = 0; j < n; ++j)
        f(i,j,double(load<T>(p+ptrdiff_t(i)*rs+ptrdiff_t(j)*cs)),i*n+j);
}
template<class T, class I, class O, class F> void sparse_entries(const SparseInput& x, anndatapp::Layout layout, F f) {
    auto* p = x.values(); auto* ix = x.indices(); auto* off = x.offsets();
    if (layout == anndatapp::Layout::csr) {
        for (size_t i = 0; i < x.rows(); ++i) {
            size_t first = size_t(load<O>(off+i*sizeof(O))), last = size_t(load<O>(off+(i+1)*sizeof(O)));
            for (size_t k = first; k < last; ++k)
                f(i,size_t(load<I>(ix+k*sizeof(I))),double(load<T>(p+k*sizeof(T))),k);
        }
    } else {
        for (size_t j = 0; j < x.columns(); ++j) {
            size_t first = size_t(load<O>(off+j*sizeof(O))), last = size_t(load<O>(off+(j+1)*sizeof(O)));
            for (size_t k = first; k < last; ++k)
                f(size_t(load<I>(ix+k*sizeof(I))),j,double(load<T>(p+k*sizeof(T))),k);
        }
    }
}
template<class T, class F> void entries(const SparseInput& x, anndatapp::Layout layout, F f) {
    if (x.index_size() == 4 && x.offset_size() == 4) return sparse_entries<T,uint32_t,uint32_t>(x,layout,f);
    if (x.index_size() == 4 && x.offset_size() == 8) return sparse_entries<T,uint32_t,uint64_t>(x,layout,f);
    if (x.index_size() == 8 && x.offset_size() == 4) return sparse_entries<T,uint64_t,uint32_t>(x,layout,f);
    return sparse_entries<T,uint64_t,uint64_t>(x,layout,f);
}
QcStatus check(const DenseInput& x, anndatapp::DType type) {
    if (x.rank() != 2) return {QcCode::invalid_argument};
    if (!width(type) || x.element_size() != width(type)) return {QcCode::unsupported_type};
    if (x.extent(0) && x.extent(1) && !x.data()) return {QcCode::invalid_argument};
    return {};
}
QcStatus check(const SparseInput& x, anndatapp::DType type, anndatapp::Layout layout) {
    if (layout != anndatapp::Layout::csr && layout != anndatapp::Layout::csc) return {QcCode::invalid_argument};
    if (!width(type) || x.value_size() != width(type)) return {QcCode::unsupported_type};
    if ((x.index_size() != 4 && x.index_size() != 8) || (x.offset_size() != 4 && x.offset_size() != 8)) return {QcCode::invalid_structure};
    if (!x.offsets() || (x.nnz() && (!x.values() || !x.indices()))) return {QcCode::invalid_argument};
    return {};
}
size_t rows(const DenseInput& x) { return x.extent(0); }
size_t columns(const DenseInput& x) { return x.extent(1); }
size_t rows(const SparseInput& x) { return x.rows(); }
size_t columns(const SparseInput& x) { return x.columns(); }
template<class F> QcStatus traverse(const DenseInput& x, anndatapp::DType t, anndatapp::Layout, F f) {
    return dispatch(t,[&]<class T>() { entries<T>(x,f); return QcStatus{}; });
}
template<class F> QcStatus traverse(const SparseInput& x, anndatapp::DType t, anndatapp::Layout l, F f) {
    return dispatch(t,[&]<class T>() { entries<T>(x,l,f); return QcStatus{}; });
}
QcStatus validate(const DenseInput& x, anndatapp::DType t, anndatapp::Layout) { return check(x,t); }
QcStatus validate(const SparseInput& x, anndatapp::DType t, anndatapp::Layout l) { return check(x,t,l); }
void sift(double* p, size_t n, size_t root) {
    while (root < n/2) {
        size_t child = 2*root+1;
        if (child+1 < n && p[child] < p[child+1]) ++child;
        if (!(p[root] < p[child])) return;
        double t = p[root]; p[root] = p[child]; p[child] = t; root = child;
    }
}
void sort(double* p, size_t n) {
    for (size_t k = n/2; k; --k) sift(p,n,k-1);
    for (size_t k = n; k > 1; --k) {
        double t = p[0]; p[0] = p[k-1]; p[k-1] = t; sift(p,k-1,0);
    }
}
template<bool Checked, class X> QcStatus factors(const X& x, anndatapp::DType t, anndatapp::Layout l, const NormalizationOptions& o, NormalizationBuffers& b) {
    auto s = validate(x,t,l); if (!s) return s;
    size_t m = rows(x), n = columns(x);
    if (b.factors_capacity < m) return {QcCode::insufficient_capacity,0,m};
    if (b.included_capacity < n) return {QcCode::insufficient_capacity,0,n};
    if (!o.target_is_set && b.scratch_capacity < m) return {QcCode::insufficient_capacity,0,m};
    if ((m && !b.factors) || (n && !b.included) || (!o.target_is_set && m && !b.scratch)) return {QcCode::invalid_argument};
    if ((o.target_is_set && (!__builtin_isfinite(o.target_sum) || o.target_sum <= 0)) || !__builtin_isfinite(o.max_fraction) || o.max_fraction < 0 || o.max_fraction > 1) return {QcCode::invalid_argument};
    for (size_t i = 0; i < m; ++i) b.factors[i] = 0;
    for (size_t j = 0; j < n; ++j) b.included[j] = 1;
    QcStatus result;
    traverse(x,t,l,[&](size_t i, size_t, double v, size_t k) {
        if constexpr (Checked) if (!__builtin_isfinite(v) || v < 0) { if (result) result = {QcCode::invalid_value,k}; return; }
        b.factors[i] += v;
    });
    if (!result) return result;
    if (o.exclude_highly_expressed) {
        traverse(x,t,l,[&](size_t i, size_t j, double v, size_t) { if (v > o.max_fraction*b.factors[i]) b.included[j] = 0; });
        for (size_t i = 0; i < m; ++i) b.factors[i] = 0;
        traverse(x,t,l,[&](size_t i, size_t j, double v, size_t) { if (b.included[j]) b.factors[i] += v; });
    }
    double target = o.target_is_set ? o.target_sum : 1;
    if (!o.target_is_set) {
        size_t count = 0;
        for (size_t i = 0; i < m; ++i) if (b.factors[i] > 0) b.scratch[count++] = b.factors[i];
        sort(b.scratch,count);
        if (count) { size_t k = count/2; target = count%2 ? b.scratch[k] : .5*b.scratch[k-1]+.5*b.scratch[k]; }
    }
    for (size_t i = 0; i < m; ++i) {
        if constexpr (Checked) if (!__builtin_isfinite(b.factors[i])) return {QcCode::invalid_value,i};
        b.factors[i] /= target;
    }
    b.target_sum = target; return {};
}
uint64_t magnitude(ptrdiff_t v) { return v < 0 ? uint64_t(-(v+1))+1 : uint64_t(v); }
uint64_t gcd(uint64_t a, uint64_t b) { while (b) { uint64_t r = a%b; a = b; b = r; } return a; }
bool writable(const DenseOutput& y) {
    size_t rows = y.extent(0), columns = y.extent(1);
    if (!rows || !columns) return true;
    uint64_t r = magnitude(y.stride_bytes(0)), c = magnitude(y.stride_bytes(1));
    if ((rows > 1 && !r) || (columns > 1 && !c)) return false;
    uint64_t g = gcd(r,c);
    return !g || !(c/g < rows && r/g < columns);
}
struct Range { uintptr_t first, last; };
template<class V> Range dense_range(const V& v) {
    size_t rows = v.extent(0), columns = v.extent(1);
    if (!rows || !columns) return {0,0};
    ptrdiff_t r = ptrdiff_t(rows-1)*v.stride_bytes(0), c = ptrdiff_t(columns-1)*v.stride_bytes(1);
    uintptr_t p = reinterpret_cast<uintptr_t>(v.data());
    return {p+(r < 0 ? r : 0)+(c < 0 ? c : 0), p+(r > 0 ? r : 0)+(c > 0 ? c : 0)+v.element_size()};
}
bool intersects(Range a, Range b) { return a.first < b.last && b.first < a.last; }
template<class P> Range bytes_range(P* p, size_t n) {
    uintptr_t a = reinterpret_cast<uintptr_t>(p); return n ? Range{a,a+n} : Range{0,0};
}
QcStatus output(const DenseInput& x, const DenseOutput& y, anndatapp::DType type, bool trusted) {
    if (y.rank() != 2 || x.extent(0) != y.extent(0) || x.extent(1) != y.extent(1)) return {QcCode::invalid_argument};
    if ((type != anndatapp::DType::float32 && type != anndatapp::DType::float64) || y.element_size() != width(type)) return {QcCode::unsupported_type};
    if (y.extent(0) && y.extent(1) && !y.data()) return {QcCode::invalid_argument};
    if (!trusted) {
        if (!writable(y)) return {QcCode::invalid_argument};
        bool exact = x.data() == y.data() && x.element_size() == y.element_size() &&
                     x.stride_bytes(0) == y.stride_bytes(0) && x.stride_bytes(1) == y.stride_bytes(1);
        if (!exact && intersects(dense_range(x),dense_range(y))) return {QcCode::invalid_argument};
    }
    return {};
}
QcStatus output(const SparseInput& x, const DenseOutput& y, anndatapp::DType type, anndatapp::Layout layout, bool trusted) {
    if (y.rank() != 2 || x.rows() != y.extent(0) || x.columns() != y.extent(1)) return {QcCode::invalid_argument};
    if ((type != anndatapp::DType::float32 && type != anndatapp::DType::float64) || y.element_size() != width(type)) return {QcCode::unsupported_type};
    if (y.extent(0) && y.extent(1) && !y.data()) return {QcCode::invalid_argument};
    if (!trusted) {
        if (!writable(y)) return {QcCode::invalid_argument};
        auto d = dense_range(y);
        size_t major = layout == anndatapp::Layout::csr ? x.rows() : x.columns();
        if (intersects(d,bytes_range(x.values(),x.nnz()*x.value_size())) ||
            intersects(d,bytes_range(x.indices(),x.nnz()*x.index_size())) ||
            intersects(d,bytes_range(x.offsets(),(major+1)*x.offset_size()))) return {QcCode::invalid_argument};
    }
    return {};
}
QcStatus output(const SparseInput& x, const SparseOutput& y, anndatapp::DType type, anndatapp::Layout layout, bool trusted = false) {
    if (x.rows() != y.rows() || x.columns() != y.columns() || x.nnz() != y.nnz() || x.index_size() != y.index_size() || x.offset_size() != y.offset_size()) return {QcCode::invalid_structure};
    if ((type != anndatapp::DType::float32 && type != anndatapp::DType::float64) || y.value_size() != width(type)) return {QcCode::unsupported_type};
    if (!y.offsets() || (y.nnz() && (!y.values() || !y.indices()))) return {QcCode::invalid_argument};
    size_t major = layout == anndatapp::Layout::csr ? x.rows() : x.columns();
    if (!trusted) {
        if (__builtin_memcmp(x.offsets(),y.offsets(),(major+1)*x.offset_size()) || (x.nnz() && __builtin_memcmp(x.indices(),y.indices(),x.nnz()*x.index_size()))) return {QcCode::invalid_structure};
        auto d = bytes_range(y.values(),y.nnz()*y.value_size());
        bool exact = x.values() == y.values() && x.value_size() == y.value_size();
        if ((!exact && intersects(d,bytes_range(x.values(),x.nnz()*x.value_size()))) ||
            intersects(d,bytes_range(x.indices(),x.nnz()*x.index_size())) ||
            intersects(d,bytes_range(x.offsets(),(major+1)*x.offset_size()))) return {QcCode::invalid_argument};
    }
    return {};
}
QcStatus output(const DenseInput& x, const DenseOutput& y, anndatapp::DType t, anndatapp::Layout, bool trusted = false) { return output(x,y,t,trusted); }

template<class T, class Op> void initialize(const DenseInput&, unsigned char*, ptrdiff_t, ptrdiff_t, Op) {}
template<class T, class Op> void initialize(const SparseInput& x, unsigned char* p, ptrdiff_t rs, ptrdiff_t cs, Op op) {
    for (size_t i = 0; i < x.rows(); ++i) for (size_t j = 0; j < x.columns(); ++j)
        store<T>(p+ptrdiff_t(i)*rs+ptrdiff_t(j)*cs,T(op(i,j,0)));
}
template<class T, class X, class Op> QcStatus write_transform(const X& x, anndatapp::DType t, anndatapp::Layout l, const DenseOutput& y, Op op) {
    auto* p = reinterpret_cast<unsigned char*>(y.data());
    ptrdiff_t rs = y.stride_bytes(0), cs = y.stride_bytes(1);
    initialize<T>(x,p,rs,cs,op);
    return traverse(x,t,l,[&](size_t i, size_t j, double v, size_t) {
        store<T>(p+ptrdiff_t(i)*rs+ptrdiff_t(j)*cs,T(op(i,j,v)));
    });
}
template<class T, class Op> QcStatus write_transform(const SparseInput& x, anndatapp::DType t, anndatapp::Layout l, const SparseOutput& y, Op op) {
    auto* p = reinterpret_cast<unsigned char*>(y.values());
    return traverse(x,t,l,[&](size_t i, size_t j, double v, size_t k) { store<T>(p+k*sizeof(T),T(op(i,j,v))); });
}
template<class X, class Y, class Op> QcStatus transform(const X& x, anndatapp::DType t, anndatapp::Layout l, const Y& y, anndatapp::DType u, Op op, bool trusted = false) {
    auto s = validate(x,t,l); if (!s) return s;
    s = output(x,y,u,l,trusted); if (!s) return s;
    if (u == anndatapp::DType::float32) return write_transform<float>(x,t,l,y,op);
    if (u == anndatapp::DType::float64) return write_transform<double>(x,t,l,y,op);
    return {QcCode::unsupported_type};
}
struct Normalize {
    const double* factors;
    double operator()(size_t i, size_t, double v) const { return factors[i] > 0 ? v/factors[i] : v; }
};
struct Log1p {
    double divisor;
    double operator()(size_t, size_t, double v) const { return __builtin_log1p(v)/divisor; }
};
struct Scale {
    const ScaleOptions& options; const double* mean; const double* stddev;
    double operator()(size_t i, size_t j, double v) const {
        if (options.observations && !options.observations[i]) return v;
        v = (v-(options.zero_center ? mean[j] : 0))/stddev[j];
        if (options.clip) {
            if (v > options.max_value) v = options.max_value;
            if (options.zero_center && v < -options.max_value) v = -options.max_value;
        }
        return v;
    }
};
template<class X, class Y> QcStatus apply_factors(const X& x, anndatapp::DType t, anndatapp::Layout l, const Y& y, anndatapp::DType u, const double* f, size_t capacity, bool trusted) {
    auto s = validate(x,t,l); if (!s) return s;
    s = output(x,y,u,l,trusted); if (!s) return s;
    size_t m = rows(x);
    if (capacity < m) return {QcCode::insufficient_capacity,0,m};
    if (m && !f) return {QcCode::invalid_argument};
    if (!trusted) for (size_t i = 0; i < m; ++i)
        if (!__builtin_isfinite(f[i]) || f[i] < 0) return {QcCode::invalid_value,i};
    return transform(x,t,l,y,u,Normalize{f},true);
}
template<class X, class Y> QcStatus normalize(const X& x, anndatapp::DType t, anndatapp::Layout l, const Y& y, anndatapp::DType u, const NormalizationOptions& o, NormalizationBuffers& b) {
    auto s = validate(x,t,l); if (!s) return s;
    s = output(x,y,u,l,o.trusted); if (!s) return s;
    s = o.trusted ? factors<false>(x,t,l,o,b) : factors<true>(x,t,l,o,b); if (!s) return s;
    return apply_factors(x,t,l,y,u,b.factors,b.factors_capacity,true);
}
template<class X, class Y> QcStatus logarithm(const X& x, anndatapp::DType t, anndatapp::Layout l, const Y& y, anndatapp::DType u, double base, bool trusted) {
    if (base != 0 && (!__builtin_isfinite(base) || base <= 0 || base == 1)) return {QcCode::invalid_argument};
    return transform(x,t,l,y,u,Log1p{base == 0 ? 1 : __builtin_log(base)},trusted);
}
template<bool Checked> QcStatus stats(const DenseInput& x, anndatapp::DType t, anndatapp::Layout, const uint8_t* mask, double* mean, double* variance, size_t capacity, uint64_t* scratch, size_t scratch_capacity, size_t* n) {
    if constexpr (Checked) return feature_stats(x,t,mask,mean,variance,capacity,scratch,scratch_capacity,n);
    else return feature_stats_unchecked(x,t,mask,mean,variance,capacity,scratch,scratch_capacity,n);
}
template<bool Checked> QcStatus stats(const SparseInput& x, anndatapp::DType t, anndatapp::Layout l, const uint8_t* mask, double* mean, double* variance, size_t capacity, uint64_t* scratch, size_t scratch_capacity, size_t* n) {
    if constexpr (Checked) return feature_stats(x,t,l,mask,mean,variance,capacity,scratch,scratch_capacity,n);
    else return feature_stats_unchecked(x,t,l,mask,mean,variance,capacity,scratch,scratch_capacity,n);
}
QcStatus scale_options(const ScaleOptions& o, size_t m, size_t n, const double* mean, const double* sd, size_t capacity) {
    if (capacity < n) return {QcCode::insufficient_capacity,0,n};
    if (n && (!mean || !sd)) return {QcCode::invalid_argument};
    if (o.clip && (!__builtin_isfinite(o.max_value) || o.max_value < 0)) return {QcCode::invalid_argument};
    if (o.trusted) return {};
    bool selected = !o.observations && m;
    for (size_t i = 0; o.observations && i < m; ++i) {
        if (o.observations[i] > 1) return {QcCode::invalid_value,i};
        selected |= o.observations[i] != 0;
    }
    if (selected) for (size_t j = 0; j < n; ++j)
        if (!__builtin_isfinite(mean[j]) || !__builtin_isfinite(sd[j]) || sd[j] <= 0) return {QcCode::invalid_value,j};
    return {};
}
template<class X, class Y> QcStatus apply(const X& x, anndatapp::DType t, anndatapp::Layout l, const Y& y, anndatapp::DType u, const ScaleOptions& o, const double* mean, const double* sd, size_t capacity) {
    auto s = validate(x,t,l); if (!s) return s;
    s = output(x,y,u,l,o.trusted); if (!s) return s;
    s = scale_options(o,rows(x),columns(x),mean,sd,capacity); if (!s) return s;
    return transform(x,t,l,y,u,Scale{o,mean,sd},true);
}
template<class X, class Y> QcStatus fit_scale(const X& x, anndatapp::DType t, anndatapp::Layout l, const Y& y, anndatapp::DType u, const ScaleOptions& o, double* mean, double* sd, size_t capacity, uint64_t* scratch, size_t scratch_capacity) {
    auto s = validate(x,t,l); if (!s) return s;
    s = output(x,y,u,l,o.trusted); if (!s) return s;
    if (o.clip && (!__builtin_isfinite(o.max_value) || o.max_value < 0)) return {QcCode::invalid_argument};
    size_t n = 0;
    s = o.trusted ? stats<false>(x,t,l,o.observations,mean,sd,capacity,scratch,scratch_capacity,&n)
                  : stats<true>(x,t,l,o.observations,mean,sd,capacity,scratch,scratch_capacity,&n);
    if (!s) return s;
    for (size_t j = 0; j < columns(x); ++j) sd[j] = sd[j] == 0 ? 1 : __builtin_sqrt(sd[j]);
    s = scale_options(o,rows(x),columns(x),mean,sd,capacity); if (!s) return s;
    return transform(x,t,l,y,u,Scale{o,mean,sd},true);
}
} // namespace

QcStatus normalization_factors(const DenseInput& x, anndatapp::DType t, const NormalizationOptions& o, NormalizationBuffers& b) { return o.trusted ? factors<false>(x,t,anndatapp::Layout::dense,o,b) : factors<true>(x,t,anndatapp::Layout::dense,o,b); }
QcStatus normalization_factors(const SparseInput& x, anndatapp::DType t, anndatapp::Layout l, const NormalizationOptions& o, NormalizationBuffers& b) { return o.trusted ? factors<false>(x,t,l,o,b) : factors<true>(x,t,l,o,b); }
QcStatus normalize_into(const DenseInput& x, anndatapp::DType t, const DenseOutput& y, anndatapp::DType u, const NormalizationOptions& o, NormalizationBuffers& b) { return normalize(x,t,anndatapp::Layout::dense,y,u,o,b); }
QcStatus normalize_into(const SparseInput& x, anndatapp::DType t, anndatapp::Layout l, const DenseOutput& y, anndatapp::DType u, const NormalizationOptions& o, NormalizationBuffers& b) { return normalize(x,t,l,y,u,o,b); }
QcStatus normalize_into(const SparseInput& x, anndatapp::DType t, anndatapp::Layout l, const SparseOutput& y, anndatapp::DType u, anndatapp::Layout destination, const NormalizationOptions& o, NormalizationBuffers& b) { if (l != destination) return {QcCode::invalid_argument}; return normalize(x,t,l,y,u,o,b); }
QcStatus apply_normalization(const DenseInput& x, anndatapp::DType t, const DenseOutput& y, anndatapp::DType u, const double* f, size_t c, bool trusted) { return apply_factors(x,t,anndatapp::Layout::dense,y,u,f,c,trusted); }
QcStatus apply_normalization(const SparseInput& x, anndatapp::DType t, anndatapp::Layout l, const DenseOutput& y, anndatapp::DType u, const double* f, size_t c, bool trusted) { return apply_factors(x,t,l,y,u,f,c,trusted); }
QcStatus apply_normalization(const SparseInput& x, anndatapp::DType t, anndatapp::Layout l, const SparseOutput& y, anndatapp::DType u, anndatapp::Layout destination, const double* f, size_t c, bool trusted) { if (l != destination) return {QcCode::invalid_argument}; return apply_factors(x,t,l,y,u,f,c,trusted); }
QcStatus log1p_into(const DenseInput& x, anndatapp::DType t, const DenseOutput& y, anndatapp::DType u, double b, bool trusted) { return logarithm(x,t,anndatapp::Layout::dense,y,u,b,trusted); }
QcStatus log1p_into(const SparseInput& x, anndatapp::DType t, anndatapp::Layout l, const DenseOutput& y, anndatapp::DType u, double b, bool trusted) { return logarithm(x,t,l,y,u,b,trusted); }
QcStatus log1p_into(const SparseInput& x, anndatapp::DType t, anndatapp::Layout l, const SparseOutput& y, anndatapp::DType u, anndatapp::Layout destination, double b, bool trusted) { if (l != destination) return {QcCode::invalid_argument}; return logarithm(x,t,l,y,u,b,trusted); }
QcStatus apply_scale(const DenseInput& x, anndatapp::DType t, const DenseOutput& y, anndatapp::DType u, const ScaleOptions& o, const double* m, const double* d, size_t c) { return apply(x,t,anndatapp::Layout::dense,y,u,o,m,d,c); }
QcStatus apply_scale(const SparseInput& x, anndatapp::DType t, anndatapp::Layout l, const DenseOutput& y, anndatapp::DType u, const ScaleOptions& o, const double* m, const double* d, size_t c) { return apply(x,t,l,y,u,o,m,d,c); }
QcStatus apply_scale(const SparseInput& x, anndatapp::DType t, anndatapp::Layout l, const SparseOutput& y, anndatapp::DType u, anndatapp::Layout destination, const ScaleOptions& o, const double* m, const double* d, size_t c) { if (l != destination || o.zero_center) return {QcCode::invalid_argument}; return apply(x,t,l,y,u,o,m,d,c); }
QcStatus scale_into(const DenseInput& x, anndatapp::DType t, const DenseOutput& y, anndatapp::DType u, const ScaleOptions& o, double* m, double* d, size_t c, uint64_t* w, size_t wc) { return fit_scale(x,t,anndatapp::Layout::dense,y,u,o,m,d,c,w,wc); }
QcStatus scale_into(const SparseInput& x, anndatapp::DType t, anndatapp::Layout l, const DenseOutput& y, anndatapp::DType u, const ScaleOptions& o, double* m, double* d, size_t c, uint64_t* w, size_t wc) { return fit_scale(x,t,l,y,u,o,m,d,c,w,wc); }
QcStatus scale_into(const SparseInput& x, anndatapp::DType t, anndatapp::Layout l, const SparseOutput& y, anndatapp::DType u, anndatapp::Layout destination, const ScaleOptions& o, double* m, double* d, size_t c, uint64_t* w, size_t wc) { if (l != destination || o.zero_center) return {QcCode::invalid_argument}; return fit_scale(x,t,l,y,u,o,m,d,c,w,wc); }
} // namespace preprocess
