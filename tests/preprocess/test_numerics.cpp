#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <new>

import anndatapp.core;
import preprocess.qc;

namespace { bool counting = false; size_t allocations = 0; }
void* operator new(std::size_t n) { if (counting) ++allocations; if (void* p = std::malloc(n ? n : 1)) return p; throw std::bad_alloc(); }
void* operator new[](std::size_t n) { return ::operator new(n); }
void operator delete(void* p) noexcept { std::free(p); }
void operator delete[](void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }
void operator delete[](void* p, std::size_t) noexcept { std::free(p); }

namespace {
bool close(double a,double b) { return std::fabs(a-b)<1e-12; }
int fail(const char* msg) { std::fprintf(stderr,"%s\n",msg); return 1; }

template<class X> int exercise(const X& x, anndatapp::DType dtype, anndatapp::Layout layout) {
    using namespace preprocess;
    double cells[3], genes[3], mean[3], dropout[3], groups[3], pct[3], top[6], logs[3], scratch[3];
    uint64_t cell_detected[3],gene_detected[3],stored[3]; size_t observations=0;
    const uint8_t mask[3]={1,0,1}; const size_t tops[2]={1,2};
    QcOptions opt{mask,1,tops,2};
    QcBuffers out{cells,cell_detected,3,genes,gene_detected,3,mean,dropout,nullptr,nullptr,nullptr,logs,groups,pct,nullptr,3,top,6};
    allocations=0; counting=true;
    QcStatus status;
    if constexpr (requires { x.rank(); }) status=qc(x,dtype,opt,out,scratch,3);
    else status=qc(x,dtype,layout,opt,out,scratch,3);
    counting=false;
    if (!status || allocations) return fail("qc status or allocation");
    if (!close(cells[0],4)||!close(cells[2],6)||!close(genes[0],3)||!close(genes[1],4)||!close(genes[2],3)) return fail("qc totals");
    if (cell_detected[0]!=2 || cell_detected[1]!=0 || gene_detected[0]!=2 || !close(groups[0],4) || !close(groups[2],2)) return fail("qc detection/groups");
    if (!close(top[0],75)||!close(top[1],100)||!close(top[4],200.0/3)||!close(top[5],100)) return fail("qc top");
    if (!std::isnan(top[2]) || !std::isnan(pct[1])) return fail("zero row");
    double var[3];
    allocations=0; counting=true;
    if constexpr (requires { x.rank(); }) status=feature_stats(x,dtype,nullptr,mean,var,3,stored,3,&observations);
    else status=feature_stats(x,dtype,layout,nullptr,mean,var,3,stored,3,&observations);
    counting=false;
    if (!status || allocations || observations!=3) return fail("stats status or allocation");
    if (!close(mean[0],1)||!close(var[0],1)||!close(mean[1],4.0/3)||!close(var[1],16.0/3)||!close(var[2],3)) return fail("stats values");
    uint8_t keep[3]; FilterBounds b; b.has_min_counts=true; b.min_counts=5;
    allocations=0; counting=true; status=filter(cells,cell_detected,3,b,keep,3); counting=false;
    if (!status || allocations || keep[0] || keep[1] || !keep[2]) return fail("filter");
    status=qc_check(3,3,opt,out,2); if (status.code!=QcCode::insufficient_capacity || status.required!=3) return fail("scratch capacity");
    auto attempt = [&](const QcBuffers& buffers,double* work,size_t work_size) {
        if constexpr (requires { x.rank(); }) return qc_accumulate(x,dtype,0,3,opt,buffers,work,work_size);
        else return qc_accumulate(x,dtype,layout,0,3,opt,buffers,work,work_size);
    };
    auto unchanged = [&] {
        for (size_t i=0;i<3;++i)
            if (cells[i]!=31 || genes[i]!=32 || groups[i]!=33 || scratch[i]!=34 || cell_detected[i]!=35 || gene_detected[i]!=36) return false;
        for (double v : top) if (v!=37) return false;
        return true;
    };
    for (size_t i=0;i<3;++i) { cells[i]=31; genes[i]=32; groups[i]=33; scratch[i]=34; cell_detected[i]=35; gene_detected[i]=36; }
    for (double& v : top) v=37;
    QcBuffers short_out=out; short_out.cells_capacity=2;
    status=attempt(short_out,scratch,3); if (status.code!=QcCode::insufficient_capacity || status.required!=3 || !unchanged()) return fail("partial cell capacity changed outputs");
    short_out=out; short_out.genes_capacity=2;
    status=attempt(short_out,scratch,3); if (status.code!=QcCode::insufficient_capacity || status.required!=3 || !unchanged()) return fail("partial gene capacity changed outputs");
    short_out=out; short_out.group_capacity=2;
    status=attempt(short_out,scratch,3); if (status.code!=QcCode::insufficient_capacity || status.required!=3 || !unchanged()) return fail("partial group capacity changed outputs");
    short_out=out; short_out.top_capacity=5;
    status=attempt(short_out,scratch,3); if (status.code!=QcCode::insufficient_capacity || status.required!=6 || !unchanged()) return fail("partial top capacity changed outputs");
    status=attempt(out,scratch,2); if (status.code!=QcCode::insufficient_capacity || status.required!=3 || !unchanged()) return fail("partial scratch capacity changed outputs");
    status=attempt(out,nullptr,3); if (status.code!=QcCode::invalid_argument || !unchanged()) return fail("partial null scratch changed outputs");
    allocations=0; counting=true;
    if constexpr (requires { x.rank(); }) status=qc_validate(x,dtype,true);
    else status=qc_validate(x,dtype,layout,true);
    if (status) {
        qc_zero(3,3,1,out);
        if constexpr (requires { x.rank(); }) status=qc_accumulate_unchecked(x,dtype,0,3,opt,out,scratch);
        else status=qc_accumulate_unchecked(x,dtype,layout,0,3,opt,out,scratch);
        if (status) qc_finish(3,3,opt,out);
    }
    if (status) {
        if constexpr (requires { x.rank(); }) status=feature_stats_unchecked(x,dtype,nullptr,mean,var,3,stored,3,&observations);
        else status=feature_stats_unchecked(x,dtype,layout,nullptr,mean,var,3,stored,3,&observations);
    }
    counting=false;
    if (!status || allocations || !close(var[2],3) || !close(top[0],75)) return fail("unchecked stages or allocation");
    double partial_total[]={1,2,3}, part2[]={4,5,6}; uint64_t partial_detected[]={1,2,3}, part2_detected[]={4,5,6};
    allocations=0; counting=true; qc_merge(3,partial_total,partial_detected,part2,part2_detected); counting=false;
    if (allocations || partial_total[2]!=9 || partial_detected[2]!=9) return fail("partial merge");
    double ma[3],m2a[3],mb[3],m2b[3],merged_variance[3]; uint64_t counta[3],countb[3];
    StatsPartial pa{0,0,ma,m2a,3},pb{0,0,mb,m2b,3};
    allocations=0; counting=true;
    if constexpr (requires { x.rank(); }) {
        status=stats_partial(x,dtype,0,1,nullptr,pa,counta,3);
        if (status) status=stats_partial(x,dtype,1,3,nullptr,pb,countb,3);
    } else {
        status=stats_partial(x,dtype,layout,0,1,nullptr,pa,counta,3);
        if (status) status=stats_partial(x,dtype,layout,1,3,nullptr,pb,countb,3);
    }
    if (status) status=stats_merge(pa,pb);
    if (status) status=stats_finish(pa,merged_variance,3);
    counting=false;
    if (!status || allocations || pa.observations!=3 || !close(pa.mean[0],1) || !close(merged_variance[0],1) || !close(merged_variance[1],16.0/3) || !close(merged_variance[2],3)) return fail("partial stats and merge");
    const uint8_t selected[]={1,0,1};
    allocations=0; counting=true;
    if constexpr (requires { x.rank(); }) status=feature_stats(x,dtype,selected,mean,var,3,stored,3,&observations);
    else status=feature_stats(x,dtype,layout,selected,mean,var,3,stored,3,&observations);
    counting=false;
    if (!status || allocations || observations!=2 || !close(mean[0],1.5) || !close(var[0],.5) || !close(mean[1],2) || !close(var[1],8)) return fail("masked stats");
    const uint8_t none[]={0,0,0};
    if constexpr (requires { x.rank(); }) {
        status=stats_partial(x,dtype,0,1,none,pa,counta,3);
        if (status) status=stats_partial(x,dtype,1,3,none,pb,countb,3);
    } else {
        status=stats_partial(x,dtype,layout,0,1,none,pa,counta,3);
        if (status) status=stats_partial(x,dtype,layout,1,3,none,pb,countb,3);
    }
    if (status) status=stats_merge(pa,pb);
    if (status) status=stats_finish(pa,merged_variance,3);
    if (!status || pa.observations || !std::isnan(pa.mean[0]) || !std::isnan(merged_variance[0])) return fail("empty partial stats finalization");
    if constexpr (requires { x.rank(); }) status=feature_stats(x,dtype,none,mean,var,3,stored,3,&observations);
    else status=feature_stats(x,dtype,layout,none,mean,var,3,stored,3,&observations);
    if (!status || observations || !std::isnan(mean[0]) || !std::isnan(var[0])) return fail("empty stats selection");
    const uint8_t invalid_mask[]={1,3,0};
    if constexpr (requires { x.rank(); }) status=feature_stats(x,dtype,invalid_mask,mean,var,3,stored,3,&observations);
    else status=feature_stats(x,dtype,layout,invalid_mask,mean,var,3,stored,3,&observations);
    if (status.code!=QcCode::invalid_value || status.index!=1) return fail("invalid mask");
    return 0;
}
}
int main() {
    auto dense=anndatapp::Array::empty<double>({3,3}); auto edit=dense.edit().dense<double>();
    const double values[]={1,0,3,0,0,0,2,4,0}; for(int i=0;i<9;++i) edit.data()[i]=values[i];
    auto v=dense.view().bytes(); if (int r=exercise(v,anndatapp::DType::float64,anndatapp::Layout::dense)) return r;
    auto csr=anndatapp::Array::csr<double,uint32_t,uint64_t>(3,3,5); auto s=csr.edit().csr<double,uint32_t,uint64_t>();
    const double vv[]={1,0,3,2,4}; const uint32_t ix[]={0,1,2,0,1}; const uint64_t off[]={0,3,3,5};
    for(int i=0;i<5;++i) { s.values()[i]=vv[i]; s.unsafe_indices()[i]=ix[i]; }
    for(int i=0;i<4;++i) s.unsafe_offsets()[i]=off[i];
    auto sv=csr.view().sparse_bytes(); if (int r=exercise(sv,anndatapp::DType::float64,anndatapp::Layout::csr)) return r;
    s.unsafe_indices()[1]=0;
    auto invalid=preprocess::qc_validate(sv,anndatapp::DType::float64,anndatapp::Layout::csr,true);
    if (invalid.code!=preprocess::QcCode::invalid_structure || invalid.index!=1) return fail("duplicate sparse index");
    s.unsafe_indices()[1]=1;
    auto csc=anndatapp::Array::csc<float,uint64_t,uint32_t>(3,3,5); auto c=csc.edit().csc<float,uint64_t,uint32_t>();
    const float cv[]={1,2,4,3,0}; const uint64_t ci[]={0,2,2,0,2}; const uint32_t co[]={0,2,3,5};
    for(int i=0;i<5;++i) { c.values()[i]=cv[i]; c.unsafe_indices()[i]=ci[i]; }
    for(int i=0;i<4;++i) c.unsafe_offsets()[i]=co[i];
    auto cvw=csc.view().sparse_bytes(); if (int r=exercise(cvw,anndatapp::DType::float32,anndatapp::Layout::csc)) return r;
    double reverse_data[]={2,4,0,0,0,0,1,0,3};
    auto reverse=anndatapp::Array::borrow_dense(reverse_data+6,{3,3},{-3,1});
    auto rv=reverse.view().bytes(); if (int r=exercise(rv,anndatapp::DType::float64,anndatapp::Layout::dense)) return r;
    return 0;
}
