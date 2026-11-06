module;
#include <stddef.h>
#include <stdint.h>
module preprocess.qc;
namespace preprocess {
namespace {
using D = anndatapp::DType;
using L = anndatapp::Layout;
template<class T> T load(const void* p) { T x; __builtin_memcpy(&x,p,sizeof(T)); return x; }
template<class F> QcStatus numeric(D d, F f) {
    switch (d) {
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
template<class T, bool Checked = true, class F> QcStatus entries(const DenseInput& x, size_t first, size_t last, F f) {
    const auto* p = x.data(); size_t m = x.extent(1); ptrdiff_t rs = x.stride_bytes(0), cs = x.stride_bytes(1);
    for (size_t i = first; i < last; ++i) for (size_t j = 0; j < m; ++j) {
        double v = double(load<T>(p+ptrdiff_t(i)*rs+ptrdiff_t(j)*cs));
        if (!f(i,j,v)) return {QcCode::invalid_value,i*m+j};
    }
    return {};
}
template<class T, class I, class O, bool Checked = true, class F> QcStatus sparse_entries(const SparseInput& x, L layout, size_t first, size_t last, F f) {
    const auto* p = x.values(); const auto* ix = x.indices(); const auto* off = x.offsets();
    size_t rows = x.rows(), cols = x.columns(), major = layout == L::csr ? rows : cols, minor = layout == L::csr ? cols : rows;
    if (!off || (x.nnz() && (!p || !ix))) return {QcCode::invalid_structure};
    size_t prior_end = size_t(load<O>(off)); if constexpr (Checked) if (prior_end != 0) return {QcCode::invalid_structure};
    for (size_t a = 0; a < major; ++a) {
        size_t end = size_t(load<O>(off+(a+1)*sizeof(O)));
        if constexpr (Checked) if (end < prior_end || end > x.nnz()) return {QcCode::invalid_structure,a};
        size_t prior = 0;
        for (size_t k = prior_end; k < end; ++k) {
            size_t b = size_t(load<I>(ix+k*sizeof(I)));
            if constexpr (Checked) if (b >= minor || (k != prior_end && b <= prior)) return {QcCode::invalid_structure,k};
            prior = b; size_t i = layout == L::csr ? a : b, j = layout == L::csr ? b : a;
            if (i < first || i >= last) continue;
            double v = double(load<T>(p+k*sizeof(T)));
            if (!f(i,j,v)) return {QcCode::invalid_value,i*cols+j};
        }
        prior_end = end;
    }
    if constexpr (Checked) return prior_end == x.nnz() ? QcStatus{} : QcStatus{QcCode::invalid_structure,major};
    return {};
}
template<class T, bool Checked = true, class F> QcStatus entries(const SparseInput& x, L layout, size_t first, size_t last, F f) {
    if (layout != L::csr && layout != L::csc) return {QcCode::invalid_argument};
    if (x.index_size() == 4 && x.offset_size() == 4) return sparse_entries<T,uint32_t,uint32_t,Checked>(x,layout,first,last,f);
    if (x.index_size() == 4 && x.offset_size() == 8) return sparse_entries<T,uint32_t,uint64_t,Checked>(x,layout,first,last,f);
    if (x.index_size() == 8 && x.offset_size() == 4) return sparse_entries<T,uint64_t,uint32_t,Checked>(x,layout,first,last,f);
    if (x.index_size() == 8 && x.offset_size() == 8) return sparse_entries<T,uint64_t,uint64_t,Checked>(x,layout,first,last,f);
    return {QcCode::invalid_structure};
}
// Requires the structure validated by the first traversal, or caller trust.
// CSR touches only this row; CSC binary-searches the sorted row indices of
// each column. No per-row scan of unrelated stored entries occurs.
template<class T> void dense_row(const DenseInput& x, size_t row, double* out) {
    const auto* p=x.data()+ptrdiff_t(row)*x.stride_bytes(0); ptrdiff_t stride=x.stride_bytes(1);
    for (size_t j=0;j<x.extent(1);++j) out[j]=double(load<T>(p+ptrdiff_t(j)*stride));
}
template<class T, class I, class O> void sparse_row(const SparseInput& x, L layout, size_t row, double* out) {
    const auto* values=x.values(); const auto* indices=x.indices(); const auto* offsets=x.offsets();
    if (layout==L::csr) {
        size_t first=size_t(load<O>(offsets+row*sizeof(O))), last=size_t(load<O>(offsets+(row+1)*sizeof(O)));
        for (size_t k=first;k<last;++k) out[size_t(load<I>(indices+k*sizeof(I)))]=double(load<T>(values+k*sizeof(T)));
    } else for (size_t j=0;j<x.columns();++j) {
        size_t lo=size_t(load<O>(offsets+j*sizeof(O))), hi=size_t(load<O>(offsets+(j+1)*sizeof(O))), end=hi;
        while (lo<hi) { size_t k=lo+(hi-lo)/2; load<I>(indices+k*sizeof(I))<row ? lo=k+1 : hi=k; }
        if (lo<end && size_t(load<I>(indices+lo*sizeof(I)))==row) out[j]=double(load<T>(values+lo*sizeof(T)));
    }
}
template<class T> void sparse_row(const SparseInput& x, L layout, size_t row, double* out) {
    if (x.index_size()==4 && x.offset_size()==4) sparse_row<T,uint32_t,uint32_t>(x,layout,row,out);
    else if (x.index_size()==4 && x.offset_size()==8) sparse_row<T,uint32_t,uint64_t>(x,layout,row,out);
    else if (x.index_size()==8 && x.offset_size()==4) sparse_row<T,uint64_t,uint32_t>(x,layout,row,out);
    else sparse_row<T,uint64_t,uint64_t>(x,layout,row,out);
}
// In-place ascending heap sort; top fractions use suffix sums.
void sort(double* a, size_t n) {
    for (size_t start = n/2; start > 0;) {
        size_t root = --start;
        for (;;) { size_t child = root*2+1; if (child >= n) break; if (child+1 < n && a[child] < a[child+1]) ++child; if (a[root] >= a[child]) break; double t = a[root]; a[root] = a[child]; a[child] = t; root = child; }
    }
    for (size_t end = n; end > 1;) {
        --end; double t = a[0]; a[0] = a[end]; a[end] = t; size_t root = 0;
        for (;;) { size_t child = root*2+1; if (child >= end) break; if (child+1 < end && a[child] < a[child+1]) ++child; if (a[root] >= a[child]) break; t = a[root]; a[root] = a[child]; a[child] = t; root = child; }
    }
}
template<bool Checked = true, class X, class Visit, class Row> QcStatus accumulate(const X& x, size_t rows, size_t cols, size_t first, size_t last, const QcOptions& o, const QcBuffers& b, double* scratch, Visit visit, Row row) {
    if (first > last || last > rows || (first < last && cols && o.tops && !scratch)) return {QcCode::invalid_argument};
    auto st = visit(x,first,last,[&](size_t i,size_t j,double v) {
        if constexpr (Checked) if (!__builtin_isfinite(v) || v < 0) return false;
        b.cell_total[i] += v; b.gene_total[j] += v;
        b.cell_detected[i] += v != 0; b.gene_detected[j] += v != 0;
        for (size_t g = 0; g < o.groups; ++g) if (o.masks[g*cols+j]) b.group_total[g*rows+i] += v;
        return true;
    });
    if (!st || !o.tops) return st;
    for (size_t i = first; i < last; ++i) {
        for (size_t j = 0; j < cols; ++j) scratch[j] = 0;
        row(x,i,scratch);
        sort(scratch,cols); double sum = 0;
        for (size_t k = 1; k <= cols; ++k) { sum += scratch[cols-k]; for (size_t q = 0; q < o.tops; ++q) if (o.percent_top[q] == k) b.top_percent[i*o.tops+q] = 100*sum/b.cell_total[i]; }
    }
    return {};
}
template<bool Checked = true, class X, class Visit> QcStatus statistics(const X& x, size_t rows, size_t cols, const uint8_t* mask, double* mean, double* variance, size_t cap, uint64_t* stored, size_t scratch_cap, size_t* observations, Visit visit) {
    if (cap < cols) return {QcCode::insufficient_capacity,0,cols};
    if (scratch_cap < cols) return {QcCode::insufficient_capacity,0,cols};
    if (!observations || (cols && (!mean || !variance || !stored))) return {QcCode::invalid_argument};
    size_t n = 0; for (size_t i = 0; i < rows; ++i) { if constexpr (Checked) if (mask && mask[i] > 1) return {QcCode::invalid_value,i}; n += !mask || mask[i]; }
    *observations = n;
    for (size_t j = 0; j < cols; ++j) { mean[j] = 0; variance[j] = 0; stored[j] = 0; }
    auto st = visit(x,0,rows,[&](size_t i,size_t j,double v) { if constexpr (Checked) if (!__builtin_isfinite(v)) return false; if (!mask || mask[i]) mean[j] += v; return true; });
    if (!st) return st;
    for (size_t j = 0; j < cols; ++j) mean[j] = n ? mean[j]/double(n) : __builtin_nan("");
    st = visit(x,0,rows,[&](size_t i,size_t j,double v) { if (!mask || mask[i]) { variance[j] += (v-mean[j])*(v-mean[j]); ++stored[j]; } return true; });
    if (!st) return st;
    for (size_t j = 0; j < cols; ++j) variance[j] = n ? (variance[j]+double(n-stored[j])*mean[j]*mean[j])/double(n > 1 ? n-1 : 1) : __builtin_nan("");
    return {};
}
template<class X, class Visit> QcStatus partial(const X& x,size_t rows,size_t cols,size_t first,size_t last,const uint8_t* mask,StatsPartial& out,uint64_t* count,size_t capacity,Visit visit) {
    if (first > last || last > rows) return {QcCode::invalid_argument};
    if (out.capacity < cols) return {QcCode::insufficient_capacity,0,cols};
    if (capacity < cols) return {QcCode::insufficient_capacity,0,cols};
    if (cols && (!out.mean || !out.m2 || !count)) return {QcCode::invalid_argument};
    size_t n=0; for (size_t i=first;i<last;++i) { if (mask && mask[i]>1) return {QcCode::invalid_value,i}; n += !mask || mask[i]; }
    out.observations=n; out.columns=cols;
    for (size_t j=0;j<cols;++j) { out.mean[j]=0; out.m2[j]=0; count[j]=0; }
    auto st=visit(x,first,last,[&](size_t i,size_t j,double v) {
        if (!__builtin_isfinite(v)) return false;
        if (!mask || mask[i]) { auto k=++count[j]; double d=v-out.mean[j]; out.mean[j]+=d/double(k); out.m2[j]+=d*(v-out.mean[j]); }
        return true;
    });
    if (!st) return st;
    for (size_t j=0;j<cols;++j) {
        size_t k=count[j]; if (!n) { out.mean[j]=0; out.m2[j]=0; continue; }
        if (k<n) { double mean=out.mean[j]; out.m2[j]+=mean*mean*(double(k)*double(n-k)/double(n)); out.mean[j]=mean*double(k)/double(n); }
    }
    return {};
}
}
QcStatus qc_check(size_t rows,size_t cols,const QcOptions& o,const QcBuffers& b,size_t scratch) {
    if (b.cells_capacity < rows) return {QcCode::insufficient_capacity,0,rows};
    if (b.genes_capacity < cols) return {QcCode::insufficient_capacity,0,cols};
    if ((rows && (!b.cell_total || !b.cell_detected)) || (cols && (!b.gene_total || !b.gene_detected))) return {QcCode::invalid_argument};
    if ((o.groups && ((cols && !o.masks) || (rows && !b.group_total))) ||
        (o.tops && (!o.percent_top || (rows && !b.top_percent)))) return {QcCode::invalid_argument};
    if (o.groups && rows > SIZE_MAX/o.groups) return {QcCode::invalid_argument};
    if (o.tops && rows > SIZE_MAX/o.tops) return {QcCode::invalid_argument};
    if (o.groups && cols > SIZE_MAX/o.groups) return {QcCode::invalid_argument};
    if (b.group_capacity < rows*o.groups) return {QcCode::insufficient_capacity,0,rows*o.groups};
    if (b.top_capacity < rows*o.tops) return {QcCode::insufficient_capacity,0,rows*o.tops};
    if (rows && o.tops && scratch < cols) return {QcCode::insufficient_capacity,0,cols};
    for (size_t k = 0; k < o.tops; ++k) { if (!o.percent_top[k] || o.percent_top[k] > cols) return {QcCode::invalid_argument,k}; for (size_t j=0;j<k;++j) if (o.percent_top[k] == o.percent_top[j]) return {QcCode::invalid_argument,k}; }
    for (size_t k = 0; k < cols*o.groups; ++k) if (o.masks[k] > 1) return {QcCode::invalid_value,k};
    return {};
}
void qc_zero(size_t rows,size_t cols,size_t groups,const QcBuffers& b) {
    for (size_t i=0;i<rows;++i) { b.cell_total[i]=0; b.cell_detected[i]=0; }
    for (size_t j=0;j<cols;++j) { b.gene_total[j]=0; b.gene_detected[j]=0; }
    for (size_t k=0;k<rows*groups;++k) b.group_total[k]=0;
}
void qc_merge(size_t cols,double* total,uint64_t* detected,const double* a,const uint64_t* b) { for (size_t j=0;j<cols;++j) { total[j]+=a[j]; detected[j]+=b[j]; } }
void qc_finish(size_t rows,size_t cols,const QcOptions& o,const QcBuffers& b) {
    for (size_t i=0;i<rows;++i) { if (b.log_cell_total) b.log_cell_total[i]=__builtin_log1p(b.cell_total[i]); if (b.log_cell_detected) b.log_cell_detected[i]=__builtin_log1p(double(b.cell_detected[i])); }
    for (size_t j=0;j<cols;++j) { double mean=rows?b.gene_total[j]/double(rows):__builtin_nan(""); if (b.mean) b.mean[j]=mean; if (b.dropout) b.dropout[j]=rows?100*(1-double(b.gene_detected[j])/double(rows)):__builtin_nan(""); if (b.log_gene_total) b.log_gene_total[j]=__builtin_log1p(b.gene_total[j]); if (b.log_gene_mean) b.log_gene_mean[j]=__builtin_log1p(mean); }
    for (size_t g=0;g<o.groups;++g) for (size_t i=0;i<rows;++i) { size_t k=g*rows+i; if (b.group_percent) b.group_percent[k]=100*b.group_total[k]/b.cell_total[i]; if (b.group_log_total) b.group_log_total[k]=__builtin_log1p(b.group_total[k]); }
}
QcStatus qc_accumulate(const DenseInput& x,D d,size_t first,size_t last,const QcOptions& o,const QcBuffers& b,double* scratch,size_t scratch_elements) {
    if (x.rank()!=2) return {QcCode::invalid_argument}; size_t rows=x.extent(0),cols=x.extent(1);
    if (first>last || last>rows) return {QcCode::invalid_argument};
    auto st=qc_check(rows,cols,o,b,scratch_elements); if (!st) return st;
    if (first<last && cols && o.tops && !scratch) return {QcCode::invalid_argument};
    return numeric(d,[&]<class T>() { if (x.element_size()!=sizeof(T)) return QcStatus{QcCode::invalid_argument}; if (rows && cols && !x.data()) return QcStatus{QcCode::invalid_argument}; return accumulate(x,rows,cols,first,last,o,b,scratch,[](const auto& a,size_t f,size_t l,auto fn){return entries<T>(a,f,l,fn);},[](const auto& a,size_t i,double* p){dense_row<T>(a,i,p);}); });
}
QcStatus qc_accumulate(const SparseInput& x,D d,L layout,size_t first,size_t last,const QcOptions& o,const QcBuffers& b,double* scratch,size_t scratch_elements) {
    if (first>last || last>x.rows()) return {QcCode::invalid_argument};
    auto st=qc_check(x.rows(),x.columns(),o,b,scratch_elements); if (!st) return st;
    if (first<last && x.columns() && o.tops && !scratch) return {QcCode::invalid_argument};
    return numeric(d,[&]<class T>() { if (x.value_size()!=sizeof(T)) return QcStatus{QcCode::invalid_argument}; return accumulate(x,x.rows(),x.columns(),first,last,o,b,scratch,[&](const auto& a,size_t f,size_t l,auto fn){return entries<T>(a,layout,f,l,fn);},[&](const auto& a,size_t i,double* p){sparse_row<T>(a,layout,i,p);}); });
}
QcStatus qc_validate(const DenseInput& x,D d,bool nonnegative) {
    if (x.rank()!=2) return {QcCode::invalid_argument};
    return numeric(d,[&]<class T>() { if (x.element_size()!=sizeof(T) || (x.extent(0) && x.extent(1) && !x.data())) return QcStatus{QcCode::invalid_argument}; return entries<T>(x,0,x.extent(0),[&](size_t,size_t,double v) {return __builtin_isfinite(v) && (!nonnegative || v>=0);}); });
}
QcStatus qc_validate(const SparseInput& x,D d,L l,bool nonnegative) {
    return numeric(d,[&]<class T>() { if (x.value_size()!=sizeof(T)) return QcStatus{QcCode::invalid_argument}; return entries<T>(x,l,0,x.rows(),[&](size_t,size_t,double v) {return __builtin_isfinite(v) && (!nonnegative || v>=0);}); });
}
QcStatus qc_accumulate_unchecked(const DenseInput& x,D d,size_t first,size_t last,const QcOptions& o,const QcBuffers& b,double* scratch) {
    if (x.rank()!=2) return {QcCode::invalid_argument}; return numeric(d,[&]<class T>() { return accumulate<false>(x,x.extent(0),x.extent(1),first,last,o,b,scratch,[](const auto& a,size_t f,size_t l,auto fn){return entries<T,false>(a,f,l,fn);},[](const auto& a,size_t i,double* p){dense_row<T>(a,i,p);}); });
}
QcStatus qc_accumulate_unchecked(const SparseInput& x,D d,L l,size_t first,size_t last,const QcOptions& o,const QcBuffers& b,double* scratch) {
    return numeric(d,[&]<class T>() { return accumulate<false>(x,x.rows(),x.columns(),first,last,o,b,scratch,[&](const auto& a,size_t f,size_t z,auto fn){return entries<T,false>(a,l,f,z,fn);},[&](const auto& a,size_t i,double* p){sparse_row<T>(a,l,i,p);}); });
}
QcStatus qc(const DenseInput& x,D d,const QcOptions& o,const QcBuffers& b,double* scratch,size_t n) {
    if (x.rank()!=2) return {QcCode::invalid_argument}; auto st=qc_check(x.extent(0),x.extent(1),o,b,n); if (!st) return st; if (x.extent(0) && x.extent(1) && o.tops && !scratch) return {QcCode::invalid_argument}; qc_zero(x.extent(0),x.extent(1),o.groups,b); st=qc_accumulate(x,d,0,x.extent(0),o,b,scratch,n); if (st) qc_finish(x.extent(0),x.extent(1),o,b); return st;
}
QcStatus qc(const SparseInput& x,D d,L l,const QcOptions& o,const QcBuffers& b,double* scratch,size_t n) {
    auto st=qc_check(x.rows(),x.columns(),o,b,n); if (!st) return st; if (x.rows() && x.columns() && o.tops && !scratch) return {QcCode::invalid_argument}; qc_zero(x.rows(),x.columns(),o.groups,b); st=qc_accumulate(x,d,l,0,x.rows(),o,b,scratch,n); if (st) qc_finish(x.rows(),x.columns(),o,b); return st;
}
QcStatus filter(const double* totals,const uint64_t* detected,size_t n,const FilterBounds& b,uint8_t* keep,size_t cap) {
    if (cap<n) return {QcCode::insufficient_capacity,0,n}; if (n && (!totals || !detected || !keep)) return {QcCode::invalid_argument};
    if ((b.has_min_counts && !__builtin_isfinite(b.min_counts)) || (b.has_max_counts && !__builtin_isfinite(b.max_counts)) || (b.has_min_counts && b.has_max_counts && b.min_counts>b.max_counts) || (b.has_min_detected && b.has_max_detected && b.min_detected>b.max_detected)) return {QcCode::invalid_argument};
    for (size_t i=0;i<n;++i) keep[i]=(!b.has_min_counts || totals[i]>=b.min_counts) && (!b.has_max_counts || totals[i]<=b.max_counts) && (!b.has_min_detected || detected[i]>=b.min_detected) && (!b.has_max_detected || detected[i]<=b.max_detected);
    return {};
}
QcStatus stats_partial(const DenseInput& x,D d,size_t first,size_t last,const uint8_t* mask,StatsPartial& out,uint64_t* count,size_t capacity) {
    if (x.rank()!=2) return {QcCode::invalid_argument};
    return numeric(d,[&]<class T>() { if (x.element_size()!=sizeof(T) || (x.extent(0) && x.extent(1) && !x.data())) return QcStatus{QcCode::invalid_argument}; return partial(x,x.extent(0),x.extent(1),first,last,mask,out,count,capacity,[](const auto& a,size_t f,size_t l,auto fn) {return entries<T>(a,f,l,fn);}); });
}
QcStatus stats_partial(const SparseInput& x,D d,L l,size_t first,size_t last,const uint8_t* mask,StatsPartial& out,uint64_t* count,size_t capacity) {
    return numeric(d,[&]<class T>() { if (x.value_size()!=sizeof(T)) return QcStatus{QcCode::invalid_argument}; return partial(x,x.rows(),x.columns(),first,last,mask,out,count,capacity,[&](const auto& a,size_t f,size_t z,auto fn) {return entries<T>(a,l,f,z,fn);}); });
}
QcStatus stats_merge(StatsPartial& a,const StatsPartial& b) {
    if (a.columns!=b.columns || a.capacity<a.columns || b.capacity<b.columns || (a.columns && (!a.mean || !a.m2 || !b.mean || !b.m2))) return {QcCode::invalid_argument};
    if (a.observations>SIZE_MAX-b.observations) return {QcCode::invalid_argument};
    size_t na=a.observations,nb=b.observations,n=na+nb;
    for (size_t j=0;j<a.columns;++j) { if (!nb) continue; if (!na) {a.mean[j]=b.mean[j]; a.m2[j]=b.m2[j]; continue;} double d=b.mean[j]-a.mean[j]; a.mean[j]+=d*double(nb)/double(n); a.m2[j]+=b.m2[j]+d*d*(double(na)*double(nb)/double(n)); }
    a.observations=n; return {};
}
QcStatus stats_finish(StatsPartial& p,double* variance,size_t capacity) {
    if (capacity<p.columns) return {QcCode::insufficient_capacity,0,p.columns};
    if (p.capacity<p.columns || (p.columns && (!variance || !p.mean || !p.m2))) return {QcCode::invalid_argument};
    for (size_t j=0;j<p.columns;++j) {
        if (!p.observations) { p.mean[j]=__builtin_nan(""); variance[j]=__builtin_nan(""); }
        else variance[j]=p.m2[j]/double(p.observations>1?p.observations-1:1);
    }
    return {};
}
QcStatus feature_stats(const DenseInput& x,D d,const uint8_t* mask,double* mean,double* variance,size_t cap,uint64_t* scratch,size_t n,size_t* obs) {
    if (x.rank()!=2) return {QcCode::invalid_argument}; return numeric(d,[&]<class T>() { if (x.element_size()!=sizeof(T)) return QcStatus{QcCode::invalid_argument}; if (x.extent(0) && x.extent(1) && !x.data()) return QcStatus{QcCode::invalid_argument}; return statistics(x,x.extent(0),x.extent(1),mask,mean,variance,cap,scratch,n,obs,[](const auto& a,size_t f,size_t l,auto fn){return entries<T>(a,f,l,fn);}); });
}
QcStatus feature_stats(const SparseInput& x,D d,L l,const uint8_t* mask,double* mean,double* variance,size_t cap,uint64_t* scratch,size_t n,size_t* obs) {
    return numeric(d,[&]<class T>() { if (x.value_size()!=sizeof(T)) return QcStatus{QcCode::invalid_argument}; return statistics(x,x.rows(),x.columns(),mask,mean,variance,cap,scratch,n,obs,[&](const auto& a,size_t f,size_t z,auto fn){return entries<T>(a,l,f,z,fn);}); });
}
QcStatus feature_stats_unchecked(const DenseInput& x,D d,const uint8_t* mask,double* mean,double* variance,size_t cap,uint64_t* scratch,size_t n,size_t* obs) {
    if (x.rank()!=2) return {QcCode::invalid_argument}; return numeric(d,[&]<class T>() { return statistics<false>(x,x.extent(0),x.extent(1),mask,mean,variance,cap,scratch,n,obs,[](const auto& a,size_t f,size_t l,auto fn){return entries<T,false>(a,f,l,fn);}); });
}
QcStatus feature_stats_unchecked(const SparseInput& x,D d,L l,const uint8_t* mask,double* mean,double* variance,size_t cap,uint64_t* scratch,size_t n,size_t* obs) {
    return numeric(d,[&]<class T>() { return statistics<false>(x,x.rows(),x.columns(),mask,mean,variance,cap,scratch,n,obs,[&](const auto& a,size_t f,size_t z,auto fn){return entries<T,false>(a,l,f,z,fn);}); });
}

}
