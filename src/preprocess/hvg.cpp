module;
#include <stddef.h>
#include <stdint.h>

module preprocess.hvg;

namespace preprocess {
namespace {

bool add_size(size_t a, size_t b, size_t& c) { return !__builtin_add_overflow(a,b,&c); }
bool mul_size(size_t a, size_t b, size_t& c) { return !__builtin_mul_overflow(a,b,&c); }
size_t max_size(size_t a, size_t b) { return a > b ? a : b; }
HvgStatus fail(HvgCode c, size_t i = 0, size_t bytes = 0) { return {c,i,bytes}; }

struct Sizes { size_t reals, indices, flags, bytes; };
HvgStatus sizes(size_t g, size_t b, const HvgOptions& o, Sizes& s) {
    if (!g || !b || !o.n_bins || (o.has_top && !o.n_top_genes) ||
        (o.method != HvgMethod::seurat && o.method != HvgMethod::cell_ranger) ||
        (o.has_log_base && (!__builtin_isfinite(o.log_base) || o.log_base <= 0 || o.log_base == 1)) ||
        (!o.has_top && (__builtin_isnan(o.min_mean) || __builtin_isnan(o.max_mean) ||
            __builtin_isnan(o.min_disp) || __builtin_isnan(o.max_disp) ||
            o.min_mean >= o.max_mean || o.min_disp >= o.max_disp))) return fail(HvgCode::invalid_argument);
    size_t bg, t, edges, real_bytes, index_bytes, n;
    if (!mul_size(b,g,bg) || !add_size(o.n_bins,1,edges) ||
        !mul_size(bg,2,s.reals) || !mul_size(g,4,t) || !add_size(s.reals,t,s.reals) ||
        !add_size(s.reals,max_size(edges,21),s.reals) || !mul_size(g,5,s.indices) ||
        !add_size(s.indices,b,s.indices) || !add_size(bg,g,s.flags) ||
        !mul_size(s.reals,sizeof(double),real_bytes) || !mul_size(s.indices,sizeof(size_t),index_bytes))
        return fail(HvgCode::overflow);
    size_t alignment = max_size(alignof(double),alignof(size_t));
    if (!add_size(real_bytes,alignment-1,n)) return fail(HvgCode::overflow);
    n -= n%alignment;
    if (!add_size(n,index_bytes,n) || !add_size(n,s.flags,s.bytes)) return fail(HvgCode::overflow);
    return {};
}
HvgStatus score_sizes(size_t g, const HvgOptions& o, Sizes& s) {
    auto status=sizes(g,1,o,s); if (status.code!=HvgCode::ok) return status;
    s.reals-=2*g; --s.indices; s.flags-=g;
    size_t real_bytes, index_bytes, n, alignment=max_size(alignof(double),alignof(size_t));
    if (!mul_size(s.reals,sizeof(double),real_bytes) || !mul_size(s.indices,sizeof(size_t),index_bytes) ||
        !add_size(real_bytes,alignment-1,n)) return fail(HvgCode::overflow);
    n-=n%alignment;
    if (!add_size(n,index_bytes,n) || !add_size(n,s.flags,s.bytes)) return fail(HvgCode::overflow);
    return {};
}

struct Work {
    double *sum, *square, *mean, *disp, *norm, *tmp, *edges;
    size_t *counts, *mc, *dc, *nc, *order, *bin;
    uint8_t *expressed, *selected;
};
Work bind(void* memory, size_t g, size_t b, const Sizes& s) {
    auto* p = static_cast<unsigned char*>(memory);
    size_t bg = b*g, alignment = max_size(alignof(double),alignof(size_t));
    auto* d = reinterpret_cast<double*>(p);
    Work w{};
    w.sum=d; d+=bg; w.square=d; d+=bg; w.mean=d; d+=g; w.disp=d; d+=g;
    w.norm=d; d+=g; w.tmp=d; d+=g; w.edges=d;
    size_t n = s.reals*sizeof(double); n = (n+alignment-1)/alignment*alignment;
    auto* k = reinterpret_cast<size_t*>(p+n);
    w.counts=k; k+=b; w.mc=k; k+=g; w.dc=k; k+=g; w.nc=k; k+=g;
    w.order=k; k+=g; w.bin=k; k+=g;
    w.expressed = reinterpret_cast<uint8_t*>(k); w.selected = w.expressed+bg;
    return w;
}
Work bind_score(void* memory, size_t g, const Sizes& s, const HvgMoments& m) {
    auto* p=static_cast<unsigned char*>(memory); auto* d=reinterpret_cast<double*>(p);
    Work w{}; w.sum=m.sum; w.square=m.square; w.expressed=m.expressed; w.counts=m.counts;
    w.mean=d; d+=g; w.disp=d; d+=g; w.norm=d; d+=g; w.tmp=d; d+=g; w.edges=d;
    size_t alignment=max_size(alignof(double),alignof(size_t)), n=s.reals*sizeof(double);
    n=(n+alignment-1)/alignment*alignment; auto* k=reinterpret_cast<size_t*>(p+n);
    w.mc=k; k+=g; w.dc=k; k+=g; w.nc=k; k+=g; w.order=k; k+=g; w.bin=k; k+=g;
    w.selected=reinterpret_cast<uint8_t*>(k); return w;
}
HvgStatus check_moments(size_t g, size_t b, const HvgMoments& m) {
    size_t n=0; if (!g || !b) return fail(HvgCode::invalid_argument);
    if (!mul_size(g,b,n)) return fail(HvgCode::overflow);
    if (!m.sum || !m.square || !m.expressed || !m.counts) return fail(HvgCode::invalid_argument);
    if (m.capacity<n || m.batch_capacity<b) return fail(HvgCode::insufficient_capacity,n);
    return {};
}

void sort_values(double* a, size_t n) {
    if (n < 2) return;
    for (size_t start=n/2; start; ) {
        size_t root=--start;
        while (root < n/2) { size_t c=2*root+1; if (c+1<n && a[c]<a[c+1]) ++c;
            if (!(a[root]<a[c])) break; double t=a[root]; a[root]=a[c]; a[c]=t; root=c; }
    }
    for (size_t end=n; --end;) {
        double t=a[0]; a[0]=a[end]; a[end]=t; size_t root=0;
        while (root < end/2) { size_t c=2*root+1; if (c+1<end && a[c]<a[c+1]) ++c;
            if (!(a[root]<a[c])) break; t=a[root]; a[root]=a[c]; a[c]=t; root=c; }
    }
}
double median(double* a, size_t n) {
    sort_values(a,n); return n%2 ? a[n/2] : .5*(a[n/2-1]+a[n/2]);
}
bool better(size_t a, size_t b, const size_t* nb, const double* norm) {
    if (nb[a] != nb[b]) return nb[a] > nb[b];
    if (__builtin_isnan(norm[a]) != __builtin_isnan(norm[b])) return !__builtin_isnan(norm[a]);
    if (!__builtin_isnan(norm[a]) && norm[a] != norm[b]) return norm[a] > norm[b];
    return a < b;
}
void sort_order(size_t* a, size_t n, const size_t* nb, const double* norm) {
    if (n < 2) return;
    for (size_t start=n/2; start;) {
        size_t root=--start;
        while (root < n/2) { size_t c=2*root+1; if (c+1<n && better(a[c],a[c+1],nb,norm)) ++c;
            if (!better(a[root],a[c],nb,norm)) break; size_t t=a[root]; a[root]=a[c]; a[c]=t; root=c; }
    }
    for (size_t end=n; --end;) {
        size_t t=a[0]; a[0]=a[end]; a[end]=t; size_t root=0;
        while (root < end/2) { size_t c=2*root+1; if (c+1<end && better(a[c],a[c+1],nb,norm)) ++c;
            if (!better(a[root],a[c],nb,norm)) break; t=a[root]; a[root]=a[c]; a[c]=t; root=c; }
    }
}

HvgStatus prepare(size_t rows, size_t g, const size_t* batch, size_t b, const HvgOptions& o,
                  HvgOutput out, void* memory, size_t bytes, Work& w) {
    Sizes s{}; auto status=sizes(g,b,o,s); if (status.code != HvgCode::ok) return status;
    if (!rows || (!o.batched && b != 1) || (o.batched && !batch) ||
        !out.mean || !out.dispersion || !out.normalized_dispersion || !out.selected ||
        (o.batched && (!out.nbatches || !out.intersection))) return fail(HvgCode::invalid_argument);
    if (out.capacity < g || bytes < s.bytes) return fail(HvgCode::insufficient_capacity,g,s.bytes);
    size_t alignment=max_size(alignof(double),alignof(size_t));
    if (!memory || uintptr_t(memory)%alignment) return fail(HvgCode::invalid_argument,0,s.bytes);
    w=bind(memory,g,b,s);
    for (size_t k=0; k<b*g; ++k) { w.sum[k]=w.square[k]=0; w.expressed[k]=0; }
    for (size_t k=0; k<b; ++k) w.counts[k]=0;
    for (size_t i=0; i<rows; ++i) { size_t k=o.batched ? batch[i] : 0;
        if (k>=b) return fail(HvgCode::invalid_input,i); ++w.counts[k]; }
    for (size_t k=0; k<b; ++k) if (!w.counts[k]) return fail(HvgCode::empty_batch,k);
    for (size_t j=0; j<g; ++j) { out.mean[j]=out.dispersion[j]=out.normalized_dispersion[j]=0;
        out.selected[j]=0; w.mc[j]=w.dc[j]=w.nc[j]=0; w.order[j]=0;
        if (o.batched) { out.nbatches[j]=0; out.intersection[j]=0; } }
    return {};
}

template<class T> T load(const unsigned char* p) { T v; __builtin_memcpy(&v,p,sizeof(T)); return v; }
template<class T, bool checked> HvgStatus add(Work w, size_t g, size_t i, size_t j, size_t b, const unsigned char* p,
                                const HvgOptions& o, double multiplier) {
    double v=double(load<T>(p));
    if constexpr (checked) if (!__builtin_isfinite(v)) return fail(HvgCode::invalid_input,i*g+j);
    size_t k=b*g+j; w.expressed[k] |= v != 0;
    if (o.method == HvgMethod::seurat) v=__builtin_expm1(v*multiplier);
    if constexpr (checked) if (!__builtin_isfinite(v) || v<0 || !__builtin_isfinite(v*v)) return fail(HvgCode::invalid_input,i*g+j);
    w.sum[k]+=v; w.square[k]+=v*v;
    if constexpr (checked) if (!__builtin_isfinite(w.sum[k]) || !__builtin_isfinite(w.square[k])) return fail(HvgCode::invalid_input,i*g+j);
    return {};
}

template<class T, bool checked> HvgStatus dense(const void* values, size_t rows, size_t g, ptrdiff_t sr, ptrdiff_t sc,
                                  const size_t* batch, const HvgOptions& o, Work w) {
    if (!values) return fail(HvgCode::invalid_argument);
    const auto* p=static_cast<const unsigned char*>(values); double multiplier=o.has_log_base ? __builtin_log(o.log_base) : 1;
    for (size_t i=0; i<rows; ++i) for (size_t j=0; j<g; ++j) {
        auto status=add<T,checked>(w,g,i,j,o.batched ? batch[i] : 0,p+ptrdiff_t(i)*sr+ptrdiff_t(j)*sc,o,multiplier);
        if (status.code != HvgCode::ok) return status;
    }
    return {};
}
template<class T, class I, class O, bool checked> HvgStatus sparse(const void* values, const void* indices, const void* offsets,
    anndatapp::Layout layout, size_t rows, size_t g, size_t nnz, const size_t* batch, const HvgOptions& o, Work w) {
    if (!offsets || (nnz && (!values || !indices))) return fail(HvgCode::invalid_argument);
    const auto* v=static_cast<const unsigned char*>(values); const auto* ix=static_cast<const unsigned char*>(indices);
    const auto* off=static_cast<const unsigned char*>(offsets); bool csr=layout==anndatapp::Layout::csr;
    size_t major=csr ? rows : g, minor=csr ? g : rows; double multiplier=o.has_log_base ? __builtin_log(o.log_base) : 1;
    if constexpr (checked) if (load<O>(off)!=0) return fail(HvgCode::invalid_input);
    for (size_t a=0; a<major; ++a) {
        size_t first=size_t(load<O>(off+a*sizeof(O))), end=size_t(load<O>(off+(a+1)*sizeof(O)));
        if constexpr (checked) if (first>end || end>nnz) return fail(HvgCode::invalid_input,a);
        size_t prev=0;
        for (size_t k=first; k<end; ++k) {
            size_t z=size_t(load<I>(ix+k*sizeof(I)));
            if constexpr (checked) if (z>=minor || (k>first && z<=prev)) return fail(HvgCode::invalid_input,k);
            prev=z; size_t i=csr ? a : z, j=csr ? z : a;
            auto status=add<T,checked>(w,g,i,j,o.batched ? batch[i] : 0,v+k*sizeof(T),o,multiplier);
            if (status.code != HvgCode::ok) return status;
        }
    }
    if constexpr (checked) if (size_t(load<O>(off+major*sizeof(O)))!=nnz) return fail(HvgCode::invalid_input,major);
    return {};
}

HvgStatus finish(size_t g, size_t batches, const HvgOptions& o, HvgOutput out, Work w) {
    const double nan=__builtin_nan(""), infinity=__builtin_inf();
    for (size_t b=0; b<batches; ++b) {
        size_t included=0, first=b*g;
        for (size_t j=0; j<g; ++j) {
            w.mean[j]=w.disp[j]=w.norm[j]=0; w.selected[j]=0;
            if (o.batched && !w.expressed[first+j]) continue;
            w.order[included++]=j;
            double mean=w.sum[first+j]/double(w.counts[b]);
            double var=w.square[first+j]/double(w.counts[b])-mean*mean;
            if (w.counts[b]>1) var*=double(w.counts[b])/double(w.counts[b]-1);
            if (mean==0) mean=1e-12;
            double disp=var/mean;
            if (o.method==HvgMethod::seurat) { disp=disp==0 ? nan : __builtin_log(disp); mean=__builtin_log1p(mean); }
            if (!__builtin_isfinite(mean)) return fail(HvgCode::invalid_input,j);
            w.mean[j]=mean; w.disp[j]=disp; w.norm[j]=nan;
        }
        if (!included) return fail(HvgCode::empty_batch,b);
        size_t edge_count=0;
        if (o.method==HvgMethod::seurat) {
            double lo=w.mean[w.order[0]], hi=lo;
            for (size_t k=1; k<included; ++k) { double v=w.mean[w.order[k]]; if (v<lo) lo=v; if (v>hi) hi=v; }
            bool constant=lo==hi;
            if (constant) { double d=lo==0 ? .001 : __builtin_fabs(lo)*.001; lo-=d; hi+=d; }
            double step=(hi-lo)/double(o.n_bins); edge_count=o.n_bins+1;
            for (size_t k=0; k<edge_count; ++k) w.edges[k]=lo+step*double(k);
            w.edges[edge_count-1]=hi; if (!constant) w.edges[0]-=(hi-lo)*.001;
        } else {
            for (size_t k=0; k<included; ++k) w.tmp[k]=w.mean[w.order[k]];
            sort_values(w.tmp,included); w.edges[edge_count++]=-infinity;
            for (int p=10; p<=100; p+=5) { double pos=double(included-1)*double(p)/100;
                size_t k=size_t(pos); double f=pos-double(k);
                w.edges[edge_count++]=k+1<included ? w.tmp[k]*(1-f)+w.tmp[k+1]*f : w.tmp[k]; }
            w.edges[edge_count++]=infinity;
            for (size_t k=1; k<edge_count; ++k)
                if (!(w.edges[k]>w.edges[k-1])) return fail(HvgCode::duplicate_bin_edges,k);
        }
        for (size_t k=0; k<included; ++k) {
            size_t j=w.order[k], lo=1, hi=edge_count;
            while (lo<hi) { size_t m=lo+(hi-lo)/2; w.edges[m]<w.mean[j] ? lo=m+1 : hi=m; }
            if (lo==edge_count) return fail(HvgCode::invalid_input,j);
            w.bin[j]=lo-1;
        }
        for (size_t bin=0; bin+1<edge_count; ++bin) {
            size_t n=0;
            for (size_t k=0; k<included; ++k) { size_t j=w.order[k];
                if (w.bin[j]==bin && !__builtin_isnan(w.disp[j])) w.tmp[n++]=w.disp[j]; }
            if (!n) continue;
            double avg=0, dev=0;
            if (o.method==HvgMethod::seurat) {
                for (size_t k=0; k<n; ++k) avg+=w.tmp[k]; avg/=double(n);
                if (n==1) { dev=avg; avg=0; }
                else { for (size_t k=0; k<n; ++k) dev+=(w.tmp[k]-avg)*(w.tmp[k]-avg); dev=__builtin_sqrt(dev/double(n-1)); }
            } else {
                avg=median(w.tmp,n);
                for (size_t k=0; k<n; ++k) w.tmp[k]=__builtin_fabs(w.tmp[k]-avg);
                dev=median(w.tmp,n)/0.6744897501960817;
            }
            for (size_t k=0; k<included; ++k) { size_t j=w.order[k];
                if (w.bin[j]==bin) w.norm[j]=(w.disp[j]-avg)/dev; }
        }
        if (o.has_top) {
            size_t n=0;
            for (size_t k=0; k<included; ++k) { size_t j=w.order[k]; if (!__builtin_isnan(w.norm[j])) w.tmp[n++]=w.norm[j]; }
            if (!n) return fail(HvgCode::no_scores,b);
            sort_values(w.tmp,n); size_t top=o.n_top_genes<n ? o.n_top_genes : n;
            double cutoff=w.tmp[n-top];
            for (size_t k=0; k<included; ++k) { size_t j=w.order[k]; w.selected[j]=!__builtin_isnan(w.norm[j]) && w.norm[j]>=cutoff; }
        } else for (size_t k=0; k<included; ++k) { size_t j=w.order[k]; double d=__builtin_isnan(w.norm[j]) ? 0 : w.norm[j];
            w.selected[j]=w.mean[j]>o.min_mean && w.mean[j]<o.max_mean && d>o.min_disp && d<o.max_disp; }
        for (size_t j=0; j<g; ++j) {
            if (!__builtin_isnan(w.mean[j])) { out.mean[j]+=w.mean[j]; ++w.mc[j]; }
            if (!__builtin_isnan(w.disp[j])) { out.dispersion[j]+=w.disp[j]; ++w.dc[j]; }
            if (!__builtin_isnan(w.norm[j])) { out.normalized_dispersion[j]+=w.norm[j]; ++w.nc[j]; }
            if (o.batched) out.nbatches[j]+=w.selected[j];
            out.selected[j]=w.selected[j];
        }
    }
    for (size_t j=0; j<g; ++j) {
        out.mean[j]=w.mc[j] ? out.mean[j]/double(w.mc[j]) : nan;
        out.dispersion[j]=w.dc[j] ? out.dispersion[j]/double(w.dc[j]) : nan;
        out.normalized_dispersion[j]=w.nc[j] ? out.normalized_dispersion[j]/double(w.nc[j]) : nan;
        if (o.batched) out.intersection[j]=out.nbatches[j]==batches;
    }
    if (o.batched && o.has_top) {
        auto status=hvg_rank_batched(out.normalized_dispersion,out.nbatches,g,o.n_top_genes,out.selected,w.order,g);
        if (status.code!=HvgCode::ok) return status;
    } else if (o.batched) for (size_t j=0; j<g; ++j) {
        if (__builtin_isnan(out.normalized_dispersion[j])) out.normalized_dispersion[j]=0;
        out.selected[j]=out.mean[j]>o.min_mean && out.mean[j]<o.max_mean &&
            out.normalized_dispersion[j]>o.min_disp && out.normalized_dispersion[j]<o.max_disp;
    }
    return {};
}

template<class F> HvgStatus dispatch(anndatapp::DType type, F f) {
    using D=anndatapp::DType;
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
    return fail(HvgCode::invalid_argument);
}
} // namespace

HvgStatus hvg_scratch_bytes(size_t genes, size_t batches, const HvgOptions& o, size_t& bytes) {
    Sizes s{}; auto status=sizes(genes,batches,o,s); bytes=status.code==HvgCode::ok ? s.bytes : 0; return status;
}
HvgStatus hvg_score_scratch_bytes(size_t genes, const HvgOptions& o, size_t& bytes) {
    Sizes s{}; auto status=score_sizes(genes,o,s); bytes=status.code==HvgCode::ok ? s.bytes : 0; return status;
}
HvgStatus hvg_moments_clear(size_t genes, size_t batches, HvgMoments m) {
    auto status=check_moments(genes,batches,m); if (status.code!=HvgCode::ok) return status;
    for (size_t k=0; k<genes*batches; ++k) { m.sum[k]=m.square[k]=0; m.expressed[k]=0; }
    for (size_t b=0; b<batches; ++b) m.counts[b]=0;
    return {};
}
HvgStatus hvg_moments_merge(size_t genes, size_t batches, HvgMoments out, const HvgMoments& source) {
    auto status=check_moments(genes,batches,out); if (status.code!=HvgCode::ok) return status;
    status=check_moments(genes,batches,source); if (status.code!=HvgCode::ok) return status;
    for (size_t b=0; b<batches; ++b)
        if (out.counts[b]>SIZE_MAX-source.counts[b]) return fail(HvgCode::overflow,b);
    for (size_t b=0; b<batches; ++b) out.counts[b]+=source.counts[b];
    for (size_t k=0; k<genes*batches; ++k) { out.sum[k]+=source.sum[k]; out.square[k]+=source.square[k];
        out.expressed[k]|=source.expressed[k]; }
    return {};
}
HvgStatus stage_prepare(size_t rows, size_t genes, const size_t* batch, size_t batches,
                        const HvgOptions& o, HvgMoments m, Work& w) {
    Sizes s{}; auto status=sizes(genes,batches,o,s); if (status.code!=HvgCode::ok) return status;
    status=check_moments(genes,batches,m); if (status.code!=HvgCode::ok) return status;
    if ((!o.batched && batches!=1) || (o.batched && rows && !batch)) return fail(HvgCode::invalid_argument);
    for (size_t i=0; i<rows; ++i) { size_t b=o.batched ? batch[i] : 0;
        if (b>=batches) return fail(HvgCode::invalid_input,i);
        if (m.counts[b]==SIZE_MAX) return fail(HvgCode::overflow,i);
        ++m.counts[b]; }
    w.sum=m.sum; w.square=m.square; w.expressed=m.expressed; w.counts=m.counts;
    return {};
}
template<bool checked> HvgStatus hvg_accumulate_dense_impl(const void* values, anndatapp::DType type,
    size_t rows, size_t genes, ptrdiff_t sr, ptrdiff_t sc, const size_t* batch, size_t batches,
    const HvgOptions& o, HvgMoments m) {
    Work w{}; auto status=stage_prepare(rows,genes,batch,batches,o,m,w);
    if (status.code!=HvgCode::ok || !rows) return status;
    return dispatch(type,[&]<class T>() { return dense<T,checked>(values,rows,genes,sr,sc,batch,o,w); });
}
template<bool checked> HvgStatus hvg_accumulate_sparse_impl(const void* values, anndatapp::DType type,
    const void* indices, const void* offsets, size_t iw, size_t ow, anndatapp::Layout layout,
    size_t rows, size_t genes, size_t nnz, const size_t* batch, size_t batches, const HvgOptions& o, HvgMoments m) {
    if ((layout!=anndatapp::Layout::csr && layout!=anndatapp::Layout::csc) || (iw!=4 && iw!=8) || (ow!=4 && ow!=8))
        return fail(HvgCode::invalid_argument);
    Work w{}; auto status=stage_prepare(rows,genes,batch,batches,o,m,w);
    if (status.code!=HvgCode::ok || !rows) return status;
    return dispatch(type,[&]<class T>() {
        if (iw==4 && ow==4) return sparse<T,uint32_t,uint32_t,checked>(values,indices,offsets,layout,rows,genes,nnz,batch,o,w);
        if (iw==4) return sparse<T,uint32_t,uint64_t,checked>(values,indices,offsets,layout,rows,genes,nnz,batch,o,w);
        if (ow==4) return sparse<T,uint64_t,uint32_t,checked>(values,indices,offsets,layout,rows,genes,nnz,batch,o,w);
        return sparse<T,uint64_t,uint64_t,checked>(values,indices,offsets,layout,rows,genes,nnz,batch,o,w);
    });
}
HvgStatus hvg_accumulate_dense(const void* values, anndatapp::DType type, size_t rows, size_t genes,
    ptrdiff_t sr, ptrdiff_t sc, const size_t* batch, size_t batches, const HvgOptions& o, HvgMoments m) {
    return hvg_accumulate_dense_impl<true>(values,type,rows,genes,sr,sc,batch,batches,o,m);
}
HvgStatus hvg_accumulate_dense_unchecked(const void* values, anndatapp::DType type, size_t rows, size_t genes,
    ptrdiff_t sr, ptrdiff_t sc, const size_t* batch, size_t batches, const HvgOptions& o, HvgMoments m) {
    return hvg_accumulate_dense_impl<false>(values,type,rows,genes,sr,sc,batch,batches,o,m);
}
HvgStatus hvg_accumulate_sparse(const void* values, anndatapp::DType type, const void* indices, const void* offsets,
    size_t iw, size_t ow, anndatapp::Layout layout, size_t rows, size_t genes, size_t nnz,
    const size_t* batch, size_t batches, const HvgOptions& o, HvgMoments m) {
    return hvg_accumulate_sparse_impl<true>(values,type,indices,offsets,iw,ow,layout,rows,genes,nnz,batch,batches,o,m);
}
HvgStatus hvg_accumulate_sparse_unchecked(const void* values, anndatapp::DType type, const void* indices, const void* offsets,
    size_t iw, size_t ow, anndatapp::Layout layout, size_t rows, size_t genes, size_t nnz,
    const size_t* batch, size_t batches, const HvgOptions& o, HvgMoments m) {
    return hvg_accumulate_sparse_impl<false>(values,type,indices,offsets,iw,ow,layout,rows,genes,nnz,batch,batches,o,m);
}
template<bool checked> HvgStatus hvg_score_impl(const HvgMoments& moments, size_t genes, size_t batches,
    const HvgOptions& o, HvgOutput out, void* scratch, size_t scratch_bytes) {
    Sizes s{}; auto status=score_sizes(genes,o,s); if (status.code!=HvgCode::ok) return status;
    status=check_moments(genes,batches,moments); if (status.code!=HvgCode::ok) return status;
    if ((!o.batched && batches!=1) || !out.mean || !out.dispersion || !out.normalized_dispersion || !out.selected ||
        (o.batched && (!out.nbatches || !out.intersection))) return fail(HvgCode::invalid_argument);
    if (out.capacity<genes || scratch_bytes<s.bytes) return fail(HvgCode::insufficient_capacity,genes,s.bytes);
    size_t alignment=max_size(alignof(double),alignof(size_t));
    if (!scratch || uintptr_t(scratch)%alignment) return fail(HvgCode::invalid_argument,0,s.bytes);
    for (size_t b=0; b<batches; ++b) if (!moments.counts[b]) return fail(HvgCode::empty_batch,b);
    if constexpr (checked) for (size_t k=0; k<genes*batches; ++k)
        if (!__builtin_isfinite(moments.sum[k]) || moments.sum[k]<0 ||
            !__builtin_isfinite(moments.square[k]) || moments.square[k]<0 || moments.expressed[k]>1)
            return fail(HvgCode::invalid_input,k);
    Work w=bind_score(scratch,genes,s,moments);
    for (size_t j=0; j<genes; ++j) { out.mean[j]=out.dispersion[j]=out.normalized_dispersion[j]=0;
        out.selected[j]=0; w.mc[j]=w.dc[j]=w.nc[j]=0; w.order[j]=0;
        if (o.batched) { out.nbatches[j]=0; out.intersection[j]=0; } }
    return finish(genes,batches,o,out,w);
}
HvgStatus hvg_score(const HvgMoments& m, size_t genes, size_t batches, const HvgOptions& o,
    HvgOutput out, void* scratch, size_t scratch_bytes) {
    return hvg_score_impl<true>(m,genes,batches,o,out,scratch,scratch_bytes);
}
HvgStatus hvg_score_unchecked(const HvgMoments& m, size_t genes, size_t batches, const HvgOptions& o,
    HvgOutput out, void* scratch, size_t scratch_bytes) {
    return hvg_score_impl<false>(m,genes,batches,o,out,scratch,scratch_bytes);
}
HvgStatus hvg_rank_batched(const double* norm, const size_t* nb, size_t genes, size_t top,
                           uint8_t* selected, size_t* order, size_t capacity) {
    if (!genes || !top || !norm || !nb || !selected || !order) return fail(HvgCode::invalid_argument);
    if (capacity<genes) return fail(HvgCode::insufficient_capacity,genes);
    for (size_t j=0; j<genes; ++j) { order[j]=j; selected[j]=0; }
    sort_order(order,genes,nb,norm);
    size_t n=top<genes ? top : genes;
    for (size_t k=0; k<n; ++k) selected[order[k]]=1;
    return {};
}
template<bool checked> HvgStatus hvg_dense_impl(const void* values, anndatapp::DType type, size_t rows, size_t genes,
    ptrdiff_t sr, ptrdiff_t sc, const size_t* batch, size_t batches, const HvgOptions& o,
    HvgOutput out, void* scratch, size_t scratch_bytes) {
    Work w{}; auto status=prepare(rows,genes,batch,batches,o,out,scratch,scratch_bytes,w);
    if (status.code!=HvgCode::ok) return status;
    status=dispatch(type,[&]<class T>() { return dense<T,checked>(values,rows,genes,sr,sc,batch,o,w); });
    return status.code==HvgCode::ok ? finish(genes,batches,o,out,w) : status;
}
template<bool checked> HvgStatus hvg_sparse_impl(const void* values, anndatapp::DType type, const void* indices, const void* offsets,
    size_t iw, size_t ow, anndatapp::Layout layout, size_t rows, size_t genes, size_t nnz,
    const size_t* batch, size_t batches, const HvgOptions& o, HvgOutput out, void* scratch, size_t scratch_bytes) {
    if ((layout!=anndatapp::Layout::csr && layout!=anndatapp::Layout::csc) || (iw!=4 && iw!=8) || (ow!=4 && ow!=8))
        return fail(HvgCode::invalid_argument);
    Work w{}; auto status=prepare(rows,genes,batch,batches,o,out,scratch,scratch_bytes,w);
    if (status.code!=HvgCode::ok) return status;
    status=dispatch(type,[&]<class T>() {
        if (iw==4 && ow==4) return sparse<T,uint32_t,uint32_t,checked>(values,indices,offsets,layout,rows,genes,nnz,batch,o,w);
        if (iw==4) return sparse<T,uint32_t,uint64_t,checked>(values,indices,offsets,layout,rows,genes,nnz,batch,o,w);
        if (ow==4) return sparse<T,uint64_t,uint32_t,checked>(values,indices,offsets,layout,rows,genes,nnz,batch,o,w);
        return sparse<T,uint64_t,uint64_t,checked>(values,indices,offsets,layout,rows,genes,nnz,batch,o,w);
    });
    return status.code==HvgCode::ok ? finish(genes,batches,o,out,w) : status;
}
HvgStatus hvg_dense(const void* values, anndatapp::DType type, size_t rows, size_t genes,
    ptrdiff_t sr, ptrdiff_t sc, const size_t* batch, size_t batches, const HvgOptions& o,
    HvgOutput out, void* scratch, size_t scratch_bytes) {
    return hvg_dense_impl<true>(values,type,rows,genes,sr,sc,batch,batches,o,out,scratch,scratch_bytes);
}
HvgStatus hvg_dense_unchecked(const void* values, anndatapp::DType type, size_t rows, size_t genes,
    ptrdiff_t sr, ptrdiff_t sc, const size_t* batch, size_t batches, const HvgOptions& o,
    HvgOutput out, void* scratch, size_t scratch_bytes) {
    return hvg_dense_impl<false>(values,type,rows,genes,sr,sc,batch,batches,o,out,scratch,scratch_bytes);
}
HvgStatus hvg_sparse(const void* values, anndatapp::DType type, const void* indices, const void* offsets,
    size_t iw, size_t ow, anndatapp::Layout layout, size_t rows, size_t genes, size_t nnz,
    const size_t* batch, size_t batches, const HvgOptions& o, HvgOutput out, void* scratch, size_t scratch_bytes) {
    return hvg_sparse_impl<true>(values,type,indices,offsets,iw,ow,layout,rows,genes,nnz,batch,batches,o,out,scratch,scratch_bytes);
}
HvgStatus hvg_sparse_unchecked(const void* values, anndatapp::DType type, const void* indices, const void* offsets,
    size_t iw, size_t ow, anndatapp::Layout layout, size_t rows, size_t genes, size_t nnz,
    const size_t* batch, size_t batches, const HvgOptions& o, HvgOutput out, void* scratch, size_t scratch_bytes) {
    return hvg_sparse_impl<false>(values,type,indices,offsets,iw,ow,layout,rows,genes,nnz,batch,batches,o,out,scratch,scratch_bytes);
}
} // namespace preprocess
