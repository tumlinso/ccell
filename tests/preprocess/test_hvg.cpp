#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <new>
#include <vector>

import preprocess;

static size_t allocations = 0;
void* operator new(size_t n) { ++allocations; if (void* p=std::malloc(n ? n : 1)) return p; throw std::bad_alloc(); }
void operator delete(void* p) noexcept { std::free(p); }
void operator delete(void* p, size_t) noexcept { std::free(p); }
void* operator new[](size_t n) { ++allocations; if (void* p=std::malloc(n ? n : 1)) return p; throw std::bad_alloc(); }
void operator delete[](void* p) noexcept { std::free(p); }
void operator delete[](void* p, size_t) noexcept { std::free(p); }

static void check(bool ok, const char* msg) { if (!ok) { std::fprintf(stderr,"%s\n",msg); std::abort(); } }
static bool near(double a, double b, double tol=2e-5) { return std::abs(a-b)<tol; }
struct Result {
    std::vector<double> mean,disp,norm;
    std::vector<uint8_t> selected,intersection;
    std::vector<size_t> nbatches;
    preprocess::HvgStatus status;
    explicit Result(size_t n): mean(n),disp(n),norm(n),selected(n),intersection(n),nbatches(n) {}
};
static Result run(anndatapp::ArrayView a, const preprocess::HvgOptions& options,
                  const size_t* batch=nullptr, size_t batches=1, bool trusted=false) {
    size_t rows=a.extent(0), genes=a.extent(1), bytes=0;
    Result r(genes); size_t before=allocations;
    auto status=preprocess::hvg_scratch_bytes(genes,batches,options,bytes);
    check(allocations==before,"HVG scratch sizing allocated");
    check(status.code==preprocess::HvgCode::ok,"HVG sizing failed");
    void* scratch=::operator new(bytes);
    preprocess::HvgOutput output{r.mean.data(),r.disp.data(),r.norm.data(),r.selected.data(),r.intersection.data(),r.nbatches.data(),genes};
    auto dtype=a.dtype(); auto layout=a.layout();
    if (layout==anndatapp::Layout::dense) {
        auto x=a.bytes(); auto sr=x.stride_bytes(0), sc=x.stride_bytes(1); auto* p=x.data();
        before=allocations;
        r.status=trusted ? preprocess::hvg_dense_unchecked(p,dtype,rows,genes,sr,sc,batch,batches,options,output,scratch,bytes)
                         : preprocess::hvg_dense(p,dtype,rows,genes,sr,sc,batch,batches,options,output,scratch,bytes);
        check(before==allocations,"HVG dense allocated");
    } else {
        auto x=a.sparse_bytes(); auto* p=x.values(); auto* ix=x.indices(); auto* off=x.offsets();
        auto iw=x.index_size(), ow=x.offset_size(), nnz=x.nnz();
        before=allocations;
        r.status=trusted ? preprocess::hvg_sparse_unchecked(p,dtype,ix,off,iw,ow,layout,rows,genes,nnz,batch,batches,options,output,scratch,bytes)
                         : preprocess::hvg_sparse(p,dtype,ix,off,iw,ow,layout,rows,genes,nnz,batch,batches,options,output,scratch,bytes);
        check(before==allocations,"HVG sparse allocated");
    }
    ::operator delete(scratch); return r;
}
static void ok(const Result& r) { check(r.status.code==preprocess::HvgCode::ok,"HVG execution failed"); }
struct GeneMajorAccumulation {
    const double* x; size_t rows, genes;
    preprocess::HvgMoments moments;
    void operator()(size_t first, size_t count) {
        moments.counts[0]+=count;
        for (size_t j=0; j<genes; ++j) for (size_t i=first; i<first+count; ++i) {
            double raw=x[j*rows+i], value=std::expm1(raw);
            moments.sum[j]+=value; moments.square[j]+=value*value; moments.expressed[j]|=raw!=0;
        }
    }
};
static void accumulate_gene_major(void* context, size_t first, size_t count) {
    (*static_cast<GeneMajorAccumulation*>(context))(first,count);
}

int main() {
    const double counts[] = {0,1,2,0,4, 1,1,1,0,2, 0,2,4,3,1, 2,0,1,0,1, 0,3,1,1,2, 1,0,0,0,5};
    std::vector<double> x(30); for (int k=0; k<30; ++k) x[k]=std::log1p(counts[k]);
    auto a=anndatapp::Array::borrow_dense(x.data(),{6,5},{5,1});
    preprocess::HvgOptions o; o.n_bins=1; o.n_top_genes=2; o.has_top=true;
    auto r=run(a.view(),o); ok(r);
    const double means[]={.510826,.773190,.916291,.510826,1.252763};
    const double norm[]={-.806067,-.299972,-.049954,1.715893,-.559900};
    const bool selected[]={false,false,true,true,false};
    for (int j=0; j<5; ++j) { check(near(r.mean[j],means[j]),"single mean");
        check(near(r.norm[j],norm[j]),"single normalized dispersion"); check(bool(r.selected[j])==selected[j],"single selection"); }
    auto unchecked=run(a.view(),o,nullptr,1,true); ok(unchecked);
    for (int j=0; j<5; ++j) check(near(unchecked.norm[j],r.norm[j],1e-12),"unchecked dense");
    double sum1[5],sum2[5],sq1[5],sq2[5]; uint8_t expr1[5],expr2[5]; size_t count1[1],count2[1];
    preprocess::HvgMoments first{sum1,sq1,expr1,count1,5,1}, second{sum2,sq2,expr2,count2,5,1};
    check(preprocess::hvg_moments_clear(5,1,first).code==preprocess::HvgCode::ok,"first moments clear");
    double invalid_x[30]; for (int k=0; k<30; ++k) invalid_x[k]=x[k]; invalid_x[17]=__builtin_nan("");
    check(preprocess::hvg_accumulate_dense(invalid_x,anndatapp::DType::float64,6,5,5*sizeof(double),sizeof(double),nullptr,1,o,first).code==preprocess::HvgCode::invalid_input,"invalid partial accumulation");
    check(preprocess::hvg_moments_clear(5,1,first).code==preprocess::HvgCode::ok,"clear after failed accumulation");
    check(preprocess::hvg_moments_clear(5,1,second).code==preprocess::HvgCode::ok,"second moments clear");
    check(preprocess::hvg_accumulate_dense(x.data(),anndatapp::DType::float64,3,5,5*sizeof(double),sizeof(double),nullptr,1,o,first).code==preprocess::HvgCode::ok,"checked partial accumulation");
    check(preprocess::hvg_accumulate_dense_unchecked(x.data()+15,anndatapp::DType::float64,3,5,5*sizeof(double),sizeof(double),nullptr,1,o,second).code==preprocess::HvgCode::ok,"unchecked partial accumulation");
    check(preprocess::hvg_moments_merge(5,1,first,second).code==preprocess::HvgCode::ok,"partial merge");
    size_t score_bytes=0; check(preprocess::hvg_score_scratch_bytes(5,o,score_bytes).code==preprocess::HvgCode::ok,"score sizing");
    void* score_scratch=::operator new(score_bytes);
    double stage_mean[5],stage_disp[5],stage_norm[5]; uint8_t stage_sel[5];
    preprocess::HvgOutput stage_output{stage_mean,stage_disp,stage_norm,stage_sel,nullptr,nullptr,5};
    check(preprocess::hvg_score(first,5,1,o,stage_output,score_scratch,score_bytes).code==preprocess::HvgCode::ok,"score merged moments");
    for (int j=0; j<5; ++j) check(near(stage_norm[j],r.norm[j],1e-12) && stage_sel[j]==r.selected[j],"staged vs monolithic");
    // Replace accumulation with an entirely caller-owned feature-major layout.
    double gene_major[30]; for (int i=0; i<6; ++i) for (int j=0; j<5; ++j) gene_major[6*j+i]=x[5*i+j];
    check(preprocess::hvg_moments_clear(5,1,first).code==preprocess::HvgCode::ok,"custom moments clear");
    GeneMajorAccumulation custom{gene_major,6,5,first};
    check(preprocess::for_blocks(6,2,custom).code==preprocess::QcCode::ok,"static replacement binding");
    check(preprocess::hvg_score_unchecked(first,5,1,o,stage_output,score_scratch,score_bytes).code==preprocess::HvgCode::ok,"custom layout score");
    for (int j=0; j<5; ++j) check(near(stage_norm[j],r.norm[j],1e-12) && stage_sel[j]==r.selected[j],"static custom layout replacement");
    check(preprocess::hvg_moments_clear(5,1,first).code==preprocess::HvgCode::ok,"runtime custom moments clear");
    check(preprocess::for_blocks(6,2,accumulate_gene_major,&custom).code==preprocess::QcCode::ok,"runtime replacement binding");
    check(preprocess::hvg_score(first,5,1,o,stage_output,score_scratch,score_bytes).code==preprocess::HvgCode::ok,"runtime custom layout score");
    for (int j=0; j<5; ++j) check(near(stage_norm[j],r.norm[j],1e-12) && stage_sel[j]==r.selected[j],"runtime custom layout replacement");
    ::operator delete(score_scratch);
    size_t batches[]={0,0,0,1,1,1}; o.batched=true;
    r=run(a.view(),o,batches,2); ok(r);
    const double batch_means[]={.490415,.770223,.857399,.490415,1.251628};
    const double batch_norm[]={-.101539,-.060268,-.623373,.521206,.263974};
    const bool batch_selected[]={false,false,false,true,true};
    for (int j=0; j<5; ++j) { check(near(r.mean[j],batch_means[j]),"batch mean");
        check(near(r.norm[j],batch_norm[j]),"batch normalized dispersion"); check(bool(r.selected[j])==batch_selected[j],"batch selection"); }
    check(r.nbatches[2]==0 && !r.intersection[3],"batch diagnostics");
    double bsum[10],bsq[10],bm[5],bd[5],bn[5]; uint8_t bexpr[10],bchosen[5],bintersection[5]; size_t bcounts[2],bnsel[5];
    preprocess::HvgMoments bmoments{bsum,bsq,bexpr,bcounts,10,2};
    preprocess::HvgOutput bout{bm,bd,bn,bchosen,bintersection,bnsel,5};
    check(preprocess::hvg_moments_clear(5,2,bmoments).code==preprocess::HvgCode::ok,"batch moments clear");
    check(preprocess::hvg_accumulate_dense(x.data(),anndatapp::DType::float64,6,5,5*sizeof(double),sizeof(double),batches,2,o,bmoments).code==preprocess::HvgCode::ok,"batch accumulate");
    size_t bscore_bytes=0; check(preprocess::hvg_score_scratch_bytes(5,o,bscore_bytes).code==preprocess::HvgCode::ok,"batch score sizing");
    void* bscore_scratch=::operator new(bscore_bytes);
    check(preprocess::hvg_score(bmoments,5,2,o,bout,bscore_scratch,bscore_bytes).code==preprocess::HvgCode::ok,"batch score");
    for (int j=0; j<5; ++j) check(near(bn[j],r.norm[j],1e-12) && bchosen[j]==r.selected[j] && bnsel[j]==r.nbatches[j],"staged batched parity");
    ::operator delete(bscore_scratch);
    x[23]=0; r=run(a.view(),o,batches,2); ok(r);
    check(near(r.mean[3],.346574) && near(r.norm[3],.655359),"batch unexpressed gene");
    x[23]=std::log1p(1.0); o.batched=false;
    std::vector<double> based(30), reversed(30);
    for (int i=0; i<6; ++i) for (int j=0; j<5; ++j) { based[i*5+j]=x[i*5+j]/std::log(2.0); reversed[(5-i)*5+j]=x[i*5+j]; }
    auto base=anndatapp::Array::borrow_dense(based.data(),{6,5},{5,1}); o.has_log_base=true; o.log_base=2;
    r=run(base.view(),o); ok(r); for (int j=0; j<5; ++j) check(near(r.norm[j],norm[j]),"log base");
    o.has_log_base=false;
    auto rev=anndatapp::Array::borrow_dense(reversed.data()+25,{6,5},{-5,1});
    r=run(rev.view(),o); ok(r); for (int j=0; j<5; ++j) check(near(r.norm[j],norm[j]),"negative stride");
    std::vector<float> fp(x.begin(),x.end()); auto fa=anndatapp::Array::borrow_dense(fp.data(),{6,5},{5,1});
    r=run(fa.view(),o); ok(r); for (int j=0; j<5; ++j) check(near(r.norm[j],norm[j]),"float input");
    std::vector<double> values; std::vector<int32_t> indices; std::vector<int64_t> offsets{0};
    for (int j=0; j<5; ++j) { for (int i=0; i<6; ++i) if (x[i*5+j]!=0) { values.push_back(x[i*5+j]); indices.push_back(i); } offsets.push_back(values.size()); }
    auto csc=anndatapp::Array::borrow_sparse(anndatapp::Layout::csc,6,5,values.size(),values.data(),indices.data(),offsets.data());
    r=run(csc.view(),o); ok(r); for (int j=0; j<5; ++j) check(near(r.norm[j],norm[j]),"CSC input");
    unchecked=run(csc.view(),o,nullptr,1,true); ok(unchecked);
    for (int j=0; j<5; ++j) check(near(unchecked.norm[j],r.norm[j],1e-12),"unchecked CSC");
    check(preprocess::hvg_moments_clear(5,1,first).code==preprocess::HvgCode::ok,"sparse moments clear");
    check(preprocess::hvg_accumulate_sparse(values.data(),anndatapp::DType::float64,indices.data(),offsets.data(),4,8,
        anndatapp::Layout::csc,6,5,values.size(),nullptr,1,o,first).code==preprocess::HvgCode::ok,"sparse accumulate");
    score_scratch=::operator new(score_bytes);
    check(preprocess::hvg_score(first,5,1,o,stage_output,score_scratch,score_bytes).code==preprocess::HvgCode::ok,"sparse score");
    for (int j=0; j<5; ++j) check(near(stage_norm[j],r.norm[j],1e-12),"staged sparse parity");
    ::operator delete(score_scratch);
    values.clear(); indices.clear(); offsets.assign(1,0);
    for (int i=0; i<6; ++i) { for (int j=0; j<5; ++j) if (x[i*5+j]!=0) { values.push_back(x[i*5+j]); indices.push_back(j); } offsets.push_back(values.size()); }
    auto csr=anndatapp::Array::borrow_sparse(anndatapp::Layout::csr,6,5,values.size(),values.data(),indices.data(),offsets.data());
    r=run(csr.view(),o); ok(r); for (int j=0; j<5; ++j) check(near(r.norm[j],norm[j]),"CSR input");
    double tied[]={0,0,1, 1,1,1, 0,0,2, 2,2,2}; for (double& v:tied) v=std::log1p(v);
    auto t=anndatapp::Array::borrow_dense(tied,{4,3},{3,1}); o.n_top_genes=1;
    r=run(t.view(),o); ok(r); check(r.selected[0] && r.selected[1] && !r.selected[2],"single batch ties");
    double ties[24]; for (int k=0; k<24; ++k) ties[k]=tied[k%12];
    auto tb=anndatapp::Array::borrow_dense(ties,{8,3},{3,1}); size_t codes[]={0,0,0,0,1,1,1,1}; o.batched=true;
    r=run(tb.view(),o,codes,2); ok(r); check(r.selected[0] && !r.selected[1] && !r.selected[2],"batch tie order");
    check(r.nbatches[0]==2 && r.nbatches[1]==2,"batch tie counts"); o.batched=false;
    double ranks[]={1,1,3,__builtin_nan("")}; size_t chosen_batches[]={2,2,1,2}, order[4]; uint8_t rank_mask[4];
    size_t before_rank=allocations;
    auto rank_status=preprocess::hvg_rank_batched(ranks,chosen_batches,4,2,rank_mask,order,4);
    check(before_rank==allocations && rank_status.code==preprocess::HvgCode::ok,"batched rank allocation or status");
    check(rank_mask[0] && rank_mask[1] && !rank_mask[2] && !rank_mask[3],"batched rank order");
    rank_status=preprocess::hvg_rank_batched(ranks,chosen_batches,4,2,rank_mask,order,3);
    check(rank_status.code==preprocess::HvgCode::insufficient_capacity && rank_status.index==4,"batched rank capacity");
    std::vector<double> cr(600); for (int i=0; i<6; ++i) for (int j=0; j<100; ++j) cr[i*100+j]=(j+1)*.05+((j*17+i*11)%13)*.01;
    auto ca=anndatapp::Array::borrow_dense(cr.data(),{6,100},{100,1}); o.method=preprocess::HvgMethod::cell_ranger; o.n_top_genes=5;
    r=run(ca.view(),o); ok(r); const int probe[]={0,5,25,50,75,95};
    const double crnorm[]={4.309008168623969,-.08200719092867405,.37191876678591346,1.7138946054453674,0,1.0214361900158766};
    for (int k=0; k<6; ++k) check(near(r.norm[probe[k]],crnorm[k],1e-8),"Cell Ranger normalized dispersion");
    for (int j=0; j<100; ++j) check(bool(r.selected[j])==(j==0 || j==1 || j==40 || j==66 || j==92),"Cell Ranger selection");
    r=run(t.view(),o); check(r.status.code==preprocess::HvgCode::duplicate_bin_edges,"duplicate percentile edges");
    size_t bytes=0; o.method=preprocess::HvgMethod::seurat; check(preprocess::hvg_scratch_bytes(5,1,o,bytes).code==preprocess::HvgCode::ok,"scratch query");
    double mean[5],disp[5],norms[5]; uint8_t chosen[5]; preprocess::HvgOutput out{mean,disp,norms,chosen,nullptr,nullptr,5};
    void* scratch=::operator new(bytes);
    auto short_s=preprocess::hvg_dense(x.data(),anndatapp::DType::float64,6,5,5*sizeof(double),sizeof(double),nullptr,1,o,out,scratch,bytes-1);
    check(short_s.code==preprocess::HvgCode::insufficient_capacity && short_s.required_bytes==bytes,"scratch capacity");
    auto malformed=preprocess::hvg_dense(x.data(),anndatapp::DType::float64,6,5,5*sizeof(double),sizeof(double),nullptr,1,o,out,static_cast<char*>(scratch)+1,bytes);
    check(malformed.code==preprocess::HvgCode::invalid_argument,"scratch alignment");
    uint32_t bad_indices[]={1,0}, bad_offsets[]={0,2}; float bad_values[]={1,2};
    size_t bad_bytes=0; check(preprocess::hvg_scratch_bytes(5,1,o,bad_bytes).code==preprocess::HvgCode::ok,"sparse scratch query");
    malformed=preprocess::hvg_sparse(bad_values,anndatapp::DType::float32,bad_indices,bad_offsets,4,4,
        anndatapp::Layout::csr,1,5,2,nullptr,1,o,out,scratch,bad_bytes);
    check(malformed.code==preprocess::HvgCode::invalid_input,"unsorted sparse indices");
    size_t overflow_bytes=1;
    check(preprocess::hvg_scratch_bytes(SIZE_MAX,2,o,overflow_bytes).code==preprocess::HvgCode::overflow && !overflow_bytes,"scratch overflow");
    int32_t integer[]={0,1,2,0,4,1,1,1,0,2,0,2,4,3,1,2,0,1,0,1,0,3,1,1,2,1,0,0,0,5};
    auto ia=anndatapp::Array::borrow_dense(integer,{6,5},{5,1}); o.method=preprocess::HvgMethod::seurat; o.n_bins=1; o.n_top_genes=2;
    double integer_as_double[30]; for (int k=0; k<30; ++k) integer_as_double[k]=integer[k];
    auto da=anndatapp::Array::borrow_dense(integer_as_double,{6,5},{5,1});
    r=run(ia.view(),o); auto dr=run(da.view(),o); ok(r); ok(dr);
    for (int j=0; j<5; ++j) check(near(r.norm[j],dr.norm[j],1e-12),"integer input dispatch");
    ::operator delete(scratch); std::puts("HVG tests passed");
}
