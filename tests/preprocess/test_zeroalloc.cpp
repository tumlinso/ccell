// Allocation instrumentation is glibc-specific and observes public malloc-family
// symbols plus C++ new/new[]. It cannot see private allocators or direct mmap.
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <new>

#if defined(__has_feature)
#if __has_feature(address_sanitizer)
#define CCELL_NO_ALLOCATION_INTERPOSITION 1
#endif
#endif
#if defined(__SANITIZE_ADDRESS__) && !defined(CCELL_NO_ALLOCATION_INTERPOSITION)
#define CCELL_NO_ALLOCATION_INTERPOSITION 1
#endif

#if defined(__linux__) && defined(__GLIBC__)
import preprocess;
import anndatapp.core;

#if !defined(CCELL_NO_ALLOCATION_INTERPOSITION)
static std::atomic<bool> measuring{false};
static std::atomic<size_t> c_calls{0}, cpp_calls{0};
extern "C" void* __libc_malloc(size_t);
extern "C" void* __libc_calloc(size_t, size_t);
extern "C" void* __libc_realloc(void*, size_t);
extern "C" void __libc_free(void*);
extern "C" void* __libc_memalign(size_t, size_t);
static void count_c() { if (measuring.load(std::memory_order_relaxed)) ++c_calls; }
static void count_cpp() { if (measuring.load(std::memory_order_relaxed)) ++cpp_calls; }
extern "C" void* malloc(size_t n) noexcept { count_c(); return __libc_malloc(n); }
extern "C" void* calloc(size_t n, size_t s) noexcept { count_c(); return __libc_calloc(n,s); }
extern "C" void* realloc(void* p, size_t n) noexcept { count_c(); return __libc_realloc(p,n); }
extern "C" void free(void* p) noexcept { __libc_free(p); }
extern "C" void* aligned_alloc(size_t a, size_t n) noexcept { count_c(); return __libc_memalign(a,n); }
extern "C" int posix_memalign(void** p, size_t a, size_t n) noexcept {
    count_c(); if (a < sizeof(void*) || (a & (a-1)) || a % sizeof(void*)) return 22;
    void* q = __libc_memalign(a,n); if (!q) return 12; *p = q; return 0;
}
void* operator new(size_t n) { count_cpp(); if (void* p = malloc(n)) return p; throw std::bad_alloc(); }
void* operator new[](size_t n) { count_cpp(); if (void* p = malloc(n)) return p; throw std::bad_alloc(); }
void* operator new(size_t n, std::align_val_t a) { count_cpp(); count_c(); if (void* p = __libc_memalign(size_t(a),n)) return p; throw std::bad_alloc(); }
void* operator new[](size_t n, std::align_val_t a) { count_cpp(); count_c(); if (void* p = __libc_memalign(size_t(a),n)) return p; throw std::bad_alloc(); }
void operator delete(void* p) noexcept { free(p); }
void operator delete[](void* p) noexcept { free(p); }
void operator delete(void* p, size_t) noexcept { free(p); }
void operator delete[](void* p, size_t) noexcept { free(p); }
void operator delete(void* p, std::align_val_t) noexcept { free(p); }
void operator delete[](void* p, std::align_val_t) noexcept { free(p); }
void operator delete(void* p, size_t, std::align_val_t) noexcept { free(p); }
void operator delete[](void* p, size_t, std::align_val_t) noexcept { free(p); }
#endif

static int failures = 0;
static void add_block(void* context, size_t first, size_t count) {
    auto* total = static_cast<size_t*>(context); *total += first + count;
}
template<class F> auto measure(const char* name, F f) {
#if defined(CCELL_NO_ALLOCATION_INTERPOSITION)
    (void)name;
    return f();
#else
    auto before_c = c_calls.load(), before_cpp = cpp_calls.load();
    measuring = true; auto result = f(); measuring = false;
    auto dc = c_calls.load()-before_c, dn = cpp_calls.load()-before_cpp;
    if (dc || dn) { std::fprintf(stderr, "%s allocated: C=%zu C++=%zu\n", name, dc, dn); ++failures; }
    return result;
#endif
}
template<class F> void measure_void(const char* name, F f) {
    measure(name, [&] { f(); return 0; });
}
static void expect(bool ok, const char* name) { if (!ok) { std::fprintf(stderr,"unexpected result: %s\n",name); ++failures; } }

int main() {
    namespace a = anndatapp; namespace p = preprocess;
#if defined(CCELL_NO_ALLOCATION_INTERPOSITION)
    std::fprintf(stderr, "Allocation measurement disabled under AddressSanitizer; exercising numerical behavior only.\n");
#else
    // Verify the instrumentation before measuring the numerical calls.
    measuring = true;
    void* (*volatile m)(size_t) = &malloc;
    void* (*volatile c)(size_t,size_t) = &calloc;
    void* (*volatile r)(void*,size_t) = &realloc;
    void* sanity_c = m(17); void* sanity_cr = c(3,8); sanity_cr = r(sanity_cr,64);
    void* sanity_new = ::operator new(19); void* sanity_array = ::operator new[](23);
    void* sanity_aligned = ::operator new(64,std::align_val_t(32));
    measuring = false;
    free(sanity_c); free(sanity_cr); ::operator delete(sanity_new); ::operator delete[](sanity_array);
    ::operator delete(sanity_aligned,std::align_val_t(32));
    expect(c_calls.load() >= 6 && cpp_calls.load() >= 3, "allocator interception sanity");
#endif

    auto x = a::Array::empty<double>({4,4}); auto y = a::Array::empty<double>({4,4});
    double* values = x.edit().dense<double>().data();
    for (size_t i=0;i<16;++i) values[i] = .1*double((i*7)%13+1);
    auto dense = x.view().bytes(); auto output = y.edit().bytes();
    auto sx = a::Array::csr<double,int32_t,int64_t>(4,4,16);
    auto sy = a::Array::csr<double,int32_t,int64_t>(4,4,16);
    auto sv = sx.edit().csr<double,int32_t,int64_t>(); auto tv = sy.edit().csr<double,int32_t,int64_t>();
    auto* si = sv.unsafe_indices(); auto* so = sv.unsafe_offsets();
    auto* ti = tv.unsafe_indices(); auto* to = tv.unsafe_offsets();
    for (size_t i=0;i<4;++i) { so[i]=to[i]=int64_t(4*i); for (size_t j=0;j<4;++j) { size_t k=4*i+j; si[k]=ti[k]=int32_t(j); sv.values()[k]=values[k]; tv.values()[k]=0; } }
    so[4]=to[4]=16;
    auto sparse = sx.view().sparse_bytes(); auto sparse_out = sy.edit().sparse_bytes();
    auto empty_x = a::Array::empty<double>({0,0}); auto empty = empty_x.view().bytes();
    auto empty_y = a::Array::empty<double>({0,0}); auto empty_out = empty_y.edit().bytes();
    constexpr auto dt = a::DType::float64; constexpr auto csr = a::Layout::csr;

    double ct[4]{}, gt[4]{}, mean[4]{}, drop[4]{}, lct[4]{}, lcd[4]{}, lgt[4]{}, lgm[4]{};
    double groups[4]{}, group_pct[4]{}, group_log[4]{}, tops[8]{}, scratch[4]{}, var[4]{};
    uint64_t cd[4]{}, gd[4]{}, counts[4]{}; uint8_t mask[4]{1,0,1,0}, keep[4]{};
    size_t topn[2]{1,3};
    p::QcOptions qo{mask,1,topn,2}; p::QcBuffers qb;
    qb.cell_total=ct; qb.cell_detected=cd; qb.cells_capacity=4;
    qb.gene_total=gt; qb.gene_detected=gd; qb.genes_capacity=4; qb.mean=mean; qb.dropout=drop;
    qb.log_cell_total=lct; qb.log_cell_detected=lcd; qb.log_gene_total=lgt; qb.log_gene_mean=lgm;
    qb.group_total=groups; qb.group_percent=group_pct; qb.group_log_total=group_log; qb.group_capacity=4;
    qb.top_percent=tops; qb.top_capacity=8;
    expect(bool(measure("qc check",[&]{return p::qc_check(4,4,qo,qb,4); })),"qc check");
    expect(bool(measure("qc validate dense",[&]{return p::qc_validate(dense,dt,true);})),"qc validate dense");
    expect(bool(measure("qc validate csr",[&]{return p::qc_validate(sparse,dt,csr,true);})),"qc validate csr");
    expect(!bool(measure("qc capacity",[&]{return p::qc_check(4,4,qo,qb,0); })),"qc capacity");
    measure_void("qc zero",[&]{p::qc_zero(4,4,1,qb);});
    ct[0] = 42; scratch[0] = 47;
    auto qc_short = measure("qc partial dense capacity",[&]{return p::qc_accumulate(dense,dt,0,4,qo,qb,scratch,0);});
    expect(qc_short.code == p::QcCode::insufficient_capacity && ct[0] == 42 && scratch[0] == 47,"qc partial dense capacity preflight");
    qc_short = measure("qc partial csr capacity",[&]{return p::qc_accumulate(sparse,dt,csr,0,4,qo,qb,scratch,0);});
    expect(qc_short.code == p::QcCode::insufficient_capacity && ct[0] == 42 && scratch[0] == 47,"qc partial csr capacity preflight");
    measure_void("qc zero after short",[&]{p::qc_zero(4,4,1,qb);});
    expect(bool(measure("qc partial dense",[&]{return p::qc_accumulate(dense,dt,0,4,qo,qb,scratch,4);})),"qc partial dense");
    measure_void("qc zero sparse",[&]{p::qc_zero(4,4,1,qb);});
    expect(bool(measure("qc partial csr",[&]{return p::qc_accumulate(sparse,dt,csr,0,4,qo,qb,scratch,4);})),"qc partial csr");
    measure_void("qc zero unchecked dense",[&]{p::qc_zero(4,4,1,qb);});
    expect(bool(measure("qc unchecked dense",[&]{return p::qc_accumulate_unchecked(dense,dt,0,4,qo,qb,scratch);})),"qc unchecked dense");
    measure_void("qc zero unchecked csr",[&]{p::qc_zero(4,4,1,qb);});
    expect(bool(measure("qc unchecked csr",[&]{return p::qc_accumulate_unchecked(sparse,dt,csr,0,4,qo,qb,scratch);})),"qc unchecked csr");
    measure_void("qc finalize",[&]{p::qc_finish(4,4,qo,qb);});
    double mt[4]{}; uint64_t md[4]{};
    measure_void("qc merge",[&]{p::qc_merge(4,mt,md,gt,gd);});
    for (int i=0;i<2;++i) {
        expect(bool(measure("qc dense",[&]{return p::qc(dense,dt,qo,qb,scratch,4);})),"qc dense");
        expect(bool(measure("qc csr",[&]{return p::qc(sparse,dt,csr,qo,qb,scratch,4);})),"qc csr");
    }
    auto qs=measure("qc invalid dtype",[&]{return p::qc(dense,static_cast<a::DType>(255),qo,qb,scratch,4);}); expect(!bool(qs),"qc invalid dtype");
    p::QcOptions empty_options; p::QcBuffers empty_outputs;
    measure("qc empty",[&]{return p::qc(empty,dt,empty_options,empty_outputs,nullptr,0);});
    p::FilterBounds fb; fb.has_min_counts=true; fb.min_counts=0;
    expect(bool(measure("filter",[&]{return p::filter(ct,cd,4,fb,keep,4);})),"filter");
    expect(!bool(measure("filter capacity",[&]{return p::filter(ct,cd,4,fb,keep,0);})),"filter capacity");
    size_t observed=0;
    expect(bool(measure("stats dense",[&]{return p::feature_stats(dense,dt,nullptr,mean,var,4,counts,4,&observed);})),"stats dense");
    expect(bool(measure("stats csr",[&]{return p::feature_stats(sparse,dt,csr,nullptr,mean,var,4,counts,4,&observed);})),"stats csr");
    expect(bool(measure("stats dense unchecked",[&]{return p::feature_stats_unchecked(dense,dt,nullptr,mean,var,4,counts,4,&observed);})),"stats dense unchecked");
    expect(bool(measure("stats csr unchecked",[&]{return p::feature_stats_unchecked(sparse,dt,csr,nullptr,mean,var,4,counts,4,&observed);})),"stats csr unchecked");
    double pm1[4]{}, pm2[4]{}, pv1[4]{}, pv2[4]{};
    p::StatsPartial partial1{0,4,pm1,pv1,4}, partial2{0,4,pm2,pv2,4};
    expect(bool(measure("stats partial dense",[&]{return p::stats_partial(dense,dt,0,2,nullptr,partial1,counts,4);})),"stats partial dense");
    expect(bool(measure("stats partial csr",[&]{return p::stats_partial(sparse,dt,csr,2,4,nullptr,partial2,counts,4);})),"stats partial csr");
    expect(bool(measure("stats merge",[&]{return p::stats_merge(partial1,partial2);})),"stats merge");
    expect(bool(measure("stats finish",[&]{return p::stats_finish(partial1,var,4);})),"stats finish");
    expect(!bool(measure("stats capacity",[&]{return p::feature_stats(dense,dt,nullptr,mean,var,0,counts,4,&observed);})),"stats capacity");
    measure("stats empty",[&]{return p::feature_stats(empty,dt,nullptr,nullptr,nullptr,0,nullptr,0,&observed);});

    double factors[4]{}, nwork[4]{}, scales[4]{}; uint8_t included[4]{};
    p::NormalizationOptions no; no.target_is_set=true; no.target_sum=10;
    p::NormalizationBuffers nb{factors,4,included,4,nwork,4,1};
    expect(bool(measure("normalization prepare dense",[&]{return p::normalization_factors(dense,dt,no,nb);})),"normalization prepare dense");
    expect(bool(measure("normalization prepare csr",[&]{return p::normalization_factors(sparse,dt,csr,no,nb);})),"normalization prepare csr");
    auto automatic=no; automatic.target_is_set=false;
    expect(bool(measure("normalization automatic target",[&]{return p::normalization_factors(dense,dt,automatic,nb);})),"normalization automatic target");
    auto short_nb=nb; short_nb.scratch_capacity=0;
    expect(!bool(measure("normalization capacity",[&]{return p::normalization_factors(dense,dt,automatic,short_nb);})),"normalization capacity");
    p::NormalizationBuffers empty_nb;
    measure("normalization prepare empty",[&]{return p::normalization_factors(empty,dt,automatic,empty_nb);});
    measure("normalize empty",[&]{return p::normalize_into(empty,dt,empty_out,dt,automatic,empty_nb);});
    measure("log empty",[&]{return p::log1p_into(empty,dt,empty_out,dt);});
    for (int i=0;i<2;++i) {
        expect(bool(measure("normalize dense",[&]{return p::normalize_into(dense,dt,output,dt,no,nb);})),"normalize dense");
        expect(bool(measure("normalize csr to dense",[&]{return p::normalize_into(sparse,dt,csr,output,dt,no,nb);})),"normalize csr to dense");
        expect(bool(measure("normalize csr",[&]{return p::normalize_into(sparse,dt,csr,sparse_out,dt,csr,no,nb);})),"normalize csr");
        expect(bool(measure("log dense",[&]{return p::log1p_into(dense,dt,output,dt);})),"log dense");
        expect(bool(measure("log csr to dense",[&]{return p::log1p_into(sparse,dt,csr,output,dt);})),"log csr to dense");
        expect(bool(measure("log csr",[&]{return p::log1p_into(sparse,dt,csr,sparse_out,dt,csr);})),"log csr");
    }
    auto trusted_no=no; trusted_no.trusted=true;
    expect(bool(measure("trusted normalization prepare",[&]{return p::normalization_factors(dense,dt,trusted_no,nb);})),"trusted normalization prepare");
    expect(bool(measure("trusted normalize",[&]{return p::normalize_into(dense,dt,output,dt,trusted_no,nb);})),"trusted normalize");
    expect(bool(measure("apply normalization",[&]{return p::apply_normalization(dense,dt,output,dt,factors,4);})),"apply normalization");
    expect(bool(measure("trusted apply normalization",[&]{return p::apply_normalization(dense,dt,output,dt,factors,4,true);})),"trusted apply normalization");
    expect(bool(measure("apply normalization csr dense",[&]{return p::apply_normalization(sparse,dt,csr,output,dt,factors,4);})),"apply normalization csr dense");
    expect(bool(measure("apply normalization csr",[&]{return p::apply_normalization(sparse,dt,csr,sparse_out,dt,csr,factors,4);})),"apply normalization csr");
    expect(bool(measure("trusted log",[&]{return p::log1p_into(dense,dt,output,dt,0,true);})),"trusted log");
    double prior = tv.values()[0];
    auto wrong = measure("wrong sparse orientation",[&]{return p::normalize_into(sparse,dt,csr,sparse_out,dt,a::Layout::csc,no,nb);});
    expect(!bool(wrong) && tv.values()[0] == prior,"wrong sparse orientation preflight");
    expect(!bool(measure("log invalid base",[&]{return p::log1p_into(dense,dt,output,dt,1);})),"log invalid base");
    p::ScaleOptions sco;
    measure("scale empty",[&]{return p::scale_into(empty,dt,empty_out,dt,sco,nullptr,nullptr,0,nullptr,0);});
    expect(bool(measure("scale dense",[&]{return p::scale_into(dense,dt,output,dt,sco,mean,scales,4,counts,4);})),"scale dense");
    expect(bool(measure("scale csr dense",[&]{return p::scale_into(sparse,dt,csr,output,dt,sco,mean,scales,4,counts,4);})),"scale csr dense");
    sco.zero_center=false;
    expect(bool(measure("scale csr",[&]{return p::scale_into(sparse,dt,csr,sparse_out,dt,csr,sco,mean,scales,4,counts,4);})),"scale csr");
    expect(bool(measure("apply scale",[&]{return p::apply_scale(dense,dt,output,dt,sco,mean,scales,4);})),"apply scale");
    expect(bool(measure("apply scale csr dense",[&]{return p::apply_scale(sparse,dt,csr,output,dt,sco,mean,scales,4);})),"apply scale csr dense");
    expect(bool(measure("apply scale csr",[&]{return p::apply_scale(sparse,dt,csr,sparse_out,dt,csr,sco,mean,scales,4);})),"apply scale csr");
    auto trusted_scale=sco; trusted_scale.trusted=true;
    expect(bool(measure("trusted scale",[&]{return p::scale_into(dense,dt,output,dt,trusted_scale,mean,scales,4,counts,4);})),"trusted scale");
    expect(!bool(measure("scale capacity",[&]{return p::scale_into(dense,dt,output,dt,sco,mean,scales,0,counts,4);})),"scale capacity");

    size_t total = 0;
    auto custom = [&](size_t first, size_t count) { total += first + count; };
    expect(bool(measure("static binding",[&]{return p::for_blocks(4,2,custom);})),"static binding");
    expect(bool(measure("runtime binding",[&]{return p::for_blocks(4,2,add_block,&total);})),"runtime binding");
    expect(bool(measure("empty binding",[&]{return p::for_blocks(0,2,add_block,&total);})),"empty binding");
    expect(!bool(measure("invalid binding",[&]{return p::for_blocks(4,0,add_block,&total);})),"invalid binding");

    p::HvgOptions ho; ho.n_bins=2;
    double hm[4]{}, hd[4]{}, hn[4]{}; uint8_t hs[4]{}, hi[4]{}; size_t hb[4]{};
    p::HvgOutput hout{hm,hd,hn,hs,hi,hb,4};
    alignas(std::max_align_t) unsigned char work[65536]; size_t needed=0;
    double sum1[4]{}, square1[4]{}, sum2[4]{}, square2[4]{};
    uint8_t expressed1[4]{}, expressed2[4]{}; size_t count1[1]{}, count2[1]{};
    p::HvgMoments moments1{sum1,square1,expressed1,count1,4,1};
    p::HvgMoments moments2{sum2,square2,expressed2,count2,4,1};
    auto mh=measure("hvg moments clear",[&]{return p::hvg_moments_clear(4,1,moments1);}); expect(mh.code==p::HvgCode::ok,"hvg moments clear");
    mh=measure("hvg moments clear 2",[&]{return p::hvg_moments_clear(4,1,moments2);}); expect(mh.code==p::HvgCode::ok,"hvg moments clear 2");
    auto short_moments=moments1; short_moments.capacity=0;
    mh=measure("hvg moments capacity",[&]{return p::hvg_moments_clear(4,1,short_moments);}); expect(mh.code!=p::HvgCode::ok,"hvg moments capacity");
    mh=measure("hvg accumulate dense",[&]{return p::hvg_accumulate_dense(values,dt,4,4,4*sizeof(double),sizeof(double),nullptr,1,ho,moments1);}); expect(mh.code==p::HvgCode::ok,"hvg accumulate dense");
    mh=measure("hvg accumulate csr",[&]{return p::hvg_accumulate_sparse(sv.values(),dt,si,so,4,8,csr,4,4,16,nullptr,1,ho,moments2);}); expect(mh.code==p::HvgCode::ok,"hvg accumulate csr");
    mh=measure("hvg merge",[&]{return p::hvg_moments_merge(4,1,moments1,moments2);}); expect(mh.code==p::HvgCode::ok,"hvg merge");
    mh=measure("hvg score size",[&]{return p::hvg_score_scratch_bytes(4,ho,needed);}); expect(mh.code==p::HvgCode::ok && needed<=sizeof(work),"hvg score size");
    mh=measure("hvg score capacity",[&]{return p::hvg_score(moments1,4,1,ho,hout,work,0);}); expect(mh.code!=p::HvgCode::ok,"hvg score capacity");
    for (int i=0;i<2;++i) {
        mh=measure("hvg score",[&]{return p::hvg_score(moments1,4,1,ho,hout,work,sizeof(work));}); expect(mh.code==p::HvgCode::ok,"hvg score");
        mh=measure("hvg score unchecked",[&]{return p::hvg_score_unchecked(moments1,4,1,ho,hout,work,sizeof(work));}); expect(mh.code==p::HvgCode::ok,"hvg score unchecked");
        mh=measure("hvg accumulate dense unchecked",[&]{return p::hvg_accumulate_dense_unchecked(values,dt,4,4,4*sizeof(double),sizeof(double),nullptr,1,ho,moments1);}); expect(mh.code==p::HvgCode::ok,"hvg accumulate dense unchecked");
        mh=measure("hvg accumulate csr unchecked",[&]{return p::hvg_accumulate_sparse_unchecked(sv.values(),dt,si,so,4,8,csr,4,4,16,nullptr,1,ho,moments2);}); expect(mh.code==p::HvgCode::ok,"hvg accumulate csr unchecked");
    }
    auto hq=measure("hvg size",[&]{return p::hvg_scratch_bytes(4,1,ho,needed);});
    expect(hq.code==p::HvgCode::ok && needed<=sizeof(work),"hvg size");
    measure("hvg empty size",[&]{return p::hvg_scratch_bytes(0,0,ho,needed);});
    hq=measure("hvg size repeated",[&]{return p::hvg_scratch_bytes(4,1,ho,needed);});
    expect(hq.code==p::HvgCode::ok,"hvg size repeated");
    for (int i=0;i<2;++i) {
        auto h=measure("hvg dense",[&]{return p::hvg_dense(values,dt,4,4,4*sizeof(double),sizeof(double),nullptr,1,ho,hout,work,sizeof(work));});
        expect(h.code==p::HvgCode::ok,"hvg dense");
        h=measure("hvg csr",[&]{return p::hvg_sparse(sv.values(),dt,si,so,sizeof(int32_t),sizeof(int64_t),csr,4,4,16,nullptr,1,ho,hout,work,sizeof(work));});
        expect(h.code==p::HvgCode::ok,"hvg csr");
        h=measure("hvg dense unchecked",[&]{return p::hvg_dense_unchecked(values,dt,4,4,4*sizeof(double),sizeof(double),nullptr,1,ho,hout,work,sizeof(work));});
        expect(h.code==p::HvgCode::ok,"hvg dense unchecked");
        h=measure("hvg csr unchecked",[&]{return p::hvg_sparse_unchecked(sv.values(),dt,si,so,sizeof(int32_t),sizeof(int64_t),csr,4,4,16,nullptr,1,ho,hout,work,sizeof(work));});
        expect(h.code==p::HvgCode::ok,"hvg csr unchecked");
    }
    auto h=measure("hvg capacity",[&]{return p::hvg_dense(values,dt,4,4,4*sizeof(double),sizeof(double),nullptr,1,ho,hout,work,0);});
    expect(h.code!=p::HvgCode::ok,"hvg capacity");
    size_t batch[4]{0,1,0,1}; auto batched=ho; batched.batched=true;
    h=measure("hvg batched size",[&]{return p::hvg_scratch_bytes(4,2,batched,needed);}); expect(h.code==p::HvgCode::ok && needed<=sizeof(work),"hvg batched size");
    h=measure("hvg batched",[&]{return p::hvg_dense(values,dt,4,4,4*sizeof(double),sizeof(double),batch,2,batched,hout,work,sizeof(work));}); expect(h.code==p::HvgCode::ok,"hvg batched");
    size_t order[4]{};
    for (int i=0;i<2;++i) {
        h=measure("hvg rank",[&]{return p::hvg_rank_batched(hn,hb,4,2,hs,order,4);}); expect(h.code==p::HvgCode::ok,"hvg rank");
    }
    h=measure("hvg rank capacity",[&]{return p::hvg_rank_batched(hn,hb,4,2,hs,order,0);}); expect(h.code==p::HvgCode::insufficient_capacity,"hvg rank capacity");
    measure("hvg rank empty",[&]{return p::hvg_rank_batched(nullptr,nullptr,0,0,nullptr,nullptr,0);});
    measure("hvg sparse invalid layout",[&]{return p::hvg_sparse(sv.values(),dt,si,so,4,8,a::Layout::dense,4,4,16,nullptr,1,ho,hout,work,sizeof(work));});
    measure("hvg empty",[&]{return p::hvg_dense(nullptr,dt,0,0,0,0,nullptr,1,ho,hout,nullptr,0);});
    return failures ? 1 : 0;
}
#else
int main() { return 77; }
#endif
