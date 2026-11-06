#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <new>

import anndatapp.core;
import preprocess;

static std::size_t allocations;
void* operator new(std::size_t n) { ++allocations; if (void* p = std::malloc(n)) return p; throw std::bad_alloc(); }
void* operator new[](std::size_t n) { ++allocations; if (void* p = std::malloc(n)) return p; throw std::bad_alloc(); }
void operator delete(void* p) noexcept { std::free(p); }
void operator delete[](void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }
void operator delete[](void* p, std::size_t) noexcept { std::free(p); }

bool near(double a, double b, double e = 1e-10) { return std::abs(a-b) <= e; }

int main() {
    using D = anndatapp::DType;
    using L = anndatapp::Layout;
    auto x = anndatapp::Array::empty<double>({3,2});
    auto* xp = x.edit().dense<double>().data();
    double values[] = {1,3, 0,0, 2,6};
    for (size_t i = 0; i < 6; ++i) xp[i] = values[i];
    auto in = x.view().bytes();
    double factors[3], work[3]; uint8_t included[2];
    preprocess::NormalizationBuffers norm{factors,3,included,2,work,3,0};
    preprocess::NormalizationOptions automatic{};
    auto before = allocations;
    auto status = preprocess::normalization_factors(in,D::float64,automatic,norm);
    if (!status || allocations != before || !near(norm.target_sum,6) || !near(factors[0],2.0/3) || factors[1] != 0 || !near(factors[2],4.0/3) || !included[0] || !included[1]) return 1;
    double automatic_factors[]={factors[0],factors[1],factors[2]};
    double automatic_target=norm.target_sum;

    auto out = anndatapp::Array::empty<float>({3,2});
    auto dst = out.edit().bytes();
    preprocess::NormalizationOptions fixed{true,12,false,.05};
    before = allocations;
    status = preprocess::normalize_into(in,D::float64,dst,D::float32,fixed,norm);
    if (!status || allocations != before || !near(norm.target_sum,12) || !near(factors[0],1.0/3) || !near(factors[1],0) || !near(factors[2],2.0/3)) return 2;
    auto* yp = out.view().dense<float>().data();
    for (size_t i = 0; i < 6; ++i) if (!near(yp[i],factors[i/2] > 0 ? values[i]/factors[i/2] : values[i])) return 3;

    double tf[3], tw[3]; uint8_t ti[2];
    preprocess::NormalizationBuffers tnorm{tf,3,ti,2,tw,3,0};
    auto trusted_auto=automatic; trusted_auto.trusted=true;
    before=allocations; auto trusted_status=preprocess::normalization_factors(in,D::float64,trusted_auto,tnorm);
    if(!trusted_status || allocations!=before || !near(tnorm.target_sum,automatic_target)) return 23;
    for(size_t i=0;i<3;++i) if(!near(tf[i],automatic_factors[i])) return 24;
    auto trusted_fixed=fixed; trusted_fixed.trusted=true;
    auto out_trusted=anndatapp::Array::empty<float>({3,2}); auto dtrusted=out_trusted.edit().bytes();
    before=allocations; trusted_status=preprocess::normalize_into(in,D::float64,dtrusted,D::float32,trusted_fixed,tnorm);
    if(!trusted_status || allocations!=before || !near(tnorm.target_sum,norm.target_sum)) return 25;
    auto* ytrusted=out_trusted.view().dense<float>().data();
    for(size_t i=0;i<6;++i) if(!near(ytrusted[i],yp[i]) || !near(tf[i/2],factors[i/2])) return 26;

    double custom_factors[]={0,2,4}; auto custom_out=anndatapp::Array::empty<double>({3,2}); auto custom_dst=custom_out.edit().bytes();
    before=allocations; status=preprocess::apply_normalization(in,D::float64,custom_dst,D::float64,custom_factors,3);
    if(!status || allocations!=before) return 27;
    auto* custom=custom_out.view().dense<double>().data();
    for(size_t i=0;i<6;++i) if(!near(custom[i],custom_factors[i/2] ? values[i]/custom_factors[i/2] : values[i])) return 28;
    status=preprocess::apply_normalization(in,D::float64,custom_dst,D::float64,custom_factors,2);
    if(status.code!=preprocess::QcCode::insufficient_capacity || status.required!=3) return 29;
    double bad_factors[]={1,-1,1};
    status=preprocess::apply_normalization(in,D::float64,custom_dst,D::float64,bad_factors,3);
    if(status.code!=preprocess::QcCode::invalid_value || status.index!=1) return 30;
    preprocess::NormalizationOptions exclusion{true,1,true,.8};
    auto h = anndatapp::Array::empty<double>({1,2});
    auto* hp = h.edit().dense<double>().data(); hp[0]=99; hp[1]=1;
    auto hv = h.view().bytes(); double hf[1], hw[1]; uint8_t hi[2];
    preprocess::NormalizationBuffers hb{hf,1,hi,2,hw,1,0};
    before=allocations; status=preprocess::normalization_factors(hv,D::float64,exclusion,hb);
    if (!status || allocations != before || hi[0] || !hi[1] || !near(hf[0],1) || !near(hb.target_sum,1)) return 4;

    auto logarithm = anndatapp::Array::empty<double>({3,2}); auto ld = logarithm.edit().bytes();
    before=allocations; status=preprocess::log1p_into(in,D::float64,ld,D::float64);
    if (!status || allocations != before) return 5;
    auto* lp=logarithm.view().dense<double>().data();
    for(size_t i=0;i<6;++i) if(!near(lp[i],std::log1p(values[i]))) return 6;
    auto trusted_log=anndatapp::Array::empty<double>({3,2}); auto tld=trusted_log.edit().bytes();
    before=allocations; trusted_status=preprocess::log1p_into(in,D::float64,tld,D::float64,0,true);
    if(!trusted_status || allocations!=before) return 42;
    auto* trusted_lp=trusted_log.view().dense<double>().data();
    for(size_t i=0;i<6;++i) if(!near(trusted_lp[i],lp[i])) return 43;
    status=preprocess::log1p_into(in,D::float64,ld,D::float64,1);
    if(status.code != preprocess::QcCode::invalid_argument) return 7;

    auto scaled = anndatapp::Array::empty<double>({3,2}); auto sdout=scaled.edit().bytes();
    uint64_t stored[2]; double mean[2], deviation[2]; uint8_t selected[]={1,0,1};
    preprocess::ScaleOptions so{true,true,1.0,selected};
    before=allocations; status=preprocess::scale_into(in,D::float64,sdout,D::float64,so,mean,deviation,2,stored,2);
    if(!status || allocations!=before || !near(mean[0],1.5) || !near(mean[1],4.5) || !near(deviation[0],std::sqrt(.5)) || !near(deviation[1],std::sqrt(4.5))) return 8;
    auto* sp=scaled.view().dense<double>().data();
    if(!near(sp[0],-std::sqrt(.5)) || !near(sp[1],-std::sqrt(.5)) || sp[2]!=0 || sp[3]!=0 || !near(sp[4],std::sqrt(.5)) || !near(sp[5],std::sqrt(.5))) return 9;

    auto trusted_scale=so; trusted_scale.trusted=true;
    auto scale_trusted=anndatapp::Array::empty<double>({3,2}); auto st=scale_trusted.edit().bytes();
    double mt[2], dt[2]; uint64_t wt[2];
    before=allocations; trusted_status=preprocess::scale_into(in,D::float64,st,D::float64,trusted_scale,mt,dt,2,wt,2);
    if(!trusted_status || allocations!=before) return 31;
    auto* spt=scale_trusted.view().dense<double>().data();
    for(size_t i=0;i<6;++i) if(!near(spt[i],sp[i]) || !near(mt[i%2],mean[i%2]) || !near(dt[i%2],deviation[i%2])) return 32;

    auto rev = anndatapp::Array::borrow_dense(values+4,{3,2},{-2,1});
    auto revout = anndatapp::Array::empty<double>({3,2}); auto rd=revout.edit().bytes();
    auto ri=rev.view().bytes(); before=allocations; status=preprocess::log1p_into(ri,D::float64,rd,D::float64);
    if(!status || allocations!=before) return 10;
    auto* rp=revout.view().dense<double>().data();
    if(!near(rp[0],std::log1p(2)) || !near(rp[1],std::log1p(6)) || !near(rp[2],0) || !near(rp[3],0) || !near(rp[4],std::log1p(1)) || !near(rp[5],std::log1p(3))) return 11;

    auto csr=anndatapp::Array::csr<double,int32_t,int64_t>(3,2,4);
    auto c=csr.edit().csr<double,int32_t,int64_t>();
    double* cv=c.values(); int32_t* ci=c.unsafe_indices(); int64_t* co=c.unsafe_offsets();
    cv[0]=1;cv[1]=3;cv[2]=2;cv[3]=6;ci[0]=0;ci[1]=1;ci[2]=0;ci[3]=1;co[0]=0;co[1]=2;co[2]=2;co[3]=4;
    auto cs=csr.view().sparse_bytes();
    auto normalized_sparse=anndatapp::Array::csr<float,int32_t,int64_t>(3,2,4);
    auto nse=normalized_sparse.edit().csr<float,int32_t,int64_t>(); auto* nsi=nse.unsafe_indices(); auto* nso=nse.unsafe_offsets(); for(int i=0;i<4;++i)nsi[i]=ci[i]; for(int i=0;i<4;++i)nso[i]=co[i]; auto csout=normalized_sparse.edit().sparse_bytes();
    before=allocations; status=preprocess::normalize_into(cs,D::float64,L::csr,csout,D::float32,L::csr,fixed,norm);
    if(!status || allocations!=before || !near(factors[0],1.0/3) || !near(factors[1],0) || !near(factors[2],2.0/3)) return 12;
    auto* csv=normalized_sparse.view().csr<float,int32_t,int64_t>().values();
    if(!near(csv[0],3) || !near(csv[1],9) || !near(csv[2],3) || !near(csv[3],9)) return 13;

    auto trusted_csr=automatic; trusted_csr.target_is_set=true; trusted_csr.target_sum=12; trusted_csr.trusted=true;
    auto normalized_trusted=anndatapp::Array::csr<float,int32_t,int64_t>(3,2,4);
    auto tce=normalized_trusted.edit().csr<float,int32_t,int64_t>();
    for(int i=0;i<4;++i)tce.unsafe_indices()[i]=ci[i];
    for(int i=0;i<4;++i)tce.unsafe_offsets()[i]=co[i];
    auto tco=normalized_trusted.edit().sparse_bytes(); double cf[3],cw[3]; uint8_t csel[2];
    preprocess::NormalizationBuffers cb{cf,3,csel,2,cw,3,0};
    before=allocations; trusted_status=preprocess::normalize_into(cs,D::float64,L::csr,tco,D::float32,L::csr,trusted_csr,cb);
    if(!trusted_status || allocations!=before || !near(cb.target_sum,12)) return 33;
    auto* tcval=normalized_trusted.view().csr<float,int32_t,int64_t>().values();
    for(int i=0;i<4;++i) if(!near(tcval[i],csv[i])) return 34;
    double sparse_factors[]={0,2,4}; auto sparse_apply_out=anndatapp::Array::csr<double,int32_t,int64_t>(3,2,4);
    auto sao=sparse_apply_out.edit().csr<double,int32_t,int64_t>();
    for(int i=0;i<4;++i)sao.unsafe_indices()[i]=ci[i];
    for(int i=0;i<4;++i)sao.unsafe_offsets()[i]=co[i];
    auto saobytes=sparse_apply_out.edit().sparse_bytes(); before=allocations;
    status=preprocess::apply_normalization(cs,D::float64,L::csr,saobytes,D::float64,L::csr,sparse_factors,3);
    if(!status || allocations!=before) return 35;
    auto* sav=sparse_apply_out.view().csr<double,int32_t,int64_t>().values();
    if(!near(sav[0],1) || !near(sav[1],3) || !near(sav[2],.5) || !near(sav[3],1.5)) return 36;
    auto dense_sparse=anndatapp::Array::empty<double>({3,2}); auto dsout=dense_sparse.edit().bytes();
    before=allocations; status=preprocess::log1p_into(cs,D::float64,L::csr,dsout,D::float64);
    if(!status || allocations!=before) return 14;
    auto* dsp=dense_sparse.view().dense<double>().data();
    if(!near(dsp[0],std::log(2)) || !near(dsp[1],std::log(4)) || !near(dsp[2],0) || !near(dsp[3],0) || !near(dsp[4],std::log(3)) || !near(dsp[5],std::log(7))) return 15;
    auto trusted_sparse_log=anndatapp::Array::empty<double>({3,2}); auto tsld=trusted_sparse_log.edit().bytes();
    before=allocations; trusted_status=preprocess::log1p_into(cs,D::float64,L::csr,tsld,D::float64,0,true);
    if(!trusted_status || allocations!=before) return 44;
    auto* trusted_dsp=trusted_sparse_log.view().dense<double>().data();
    for(int i=0;i<6;++i) if(!near(trusted_dsp[i],dsp[i])) return 45;

    auto csc=anndatapp::Array::csc<double,int64_t,int32_t>(3,2,4);
    auto sc=csc.edit().csc<double,int64_t,int32_t>(); double* sv=sc.values(); int64_t* si=sc.unsafe_indices(); int32_t* sofs=sc.unsafe_offsets();
    sv[0]=1;sv[1]=2;sv[2]=3;sv[3]=6;si[0]=0;si[1]=2;si[2]=0;si[3]=2;sofs[0]=0;sofs[1]=2;sofs[2]=4;
    auto sci=csc.view().sparse_bytes();
    double csc_factors[3], csc_work[3]; uint8_t csc_included[2];
    preprocess::NormalizationBuffers csc_norm{csc_factors,3,csc_included,2,csc_work,3,0};
    auto checked_csc=automatic; before=allocations; status=preprocess::normalization_factors(sci,D::float64,L::csc,checked_csc,csc_norm);
    if(!status || allocations!=before || !near(csc_factors[0],2.0/3) || csc_factors[1]!=0 || !near(csc_factors[2],4.0/3)) return 37;
    checked_csc.trusted=true; double csc_tf[3], csc_tw[3]; uint8_t csc_ti[2];
    preprocess::NormalizationBuffers csc_tnorm{csc_tf,3,csc_ti,2,csc_tw,3,0};
    before=allocations; status=preprocess::normalization_factors(sci,D::float64,L::csc,checked_csc,csc_tnorm);
    if(!status || allocations!=before || !near(csc_tnorm.target_sum,csc_norm.target_sum)) return 38;
    for(int i=0;i<3;++i) if(!near(csc_tf[i],csc_factors[i])) return 39;
    before=allocations; status=preprocess::log1p_into(sci,D::float64,L::csc,dsout,D::float64);
    if(!status || allocations!=before) return 16;
    if(!near(dsp[0],std::log(2)) || !near(dsp[1],std::log(4)) || !near(dsp[2],0) || !near(dsp[3],0) || !near(dsp[4],std::log(3)) || !near(dsp[5],std::log(7))) return 17;
    auto trusted_csc_log=anndatapp::Array::empty<double>({3,2}); auto tccd=trusted_csc_log.edit().bytes();
    before=allocations; trusted_status=preprocess::log1p_into(sci,D::float64,L::csc,tccd,D::float64,0,true);
    if(!trusted_status || allocations!=before) return 46;
    auto* trusted_csc_dense=trusted_csc_log.view().dense<double>().data();
    for(int i=0;i<6;++i) if(!near(trusted_csc_dense[i],dsp[i])) return 47;
    auto cscout=anndatapp::Array::csc<double,int64_t,int32_t>(3,2,4);
    auto co2=cscout.edit().csc<double,int64_t,int32_t>(); for(int i=0;i<4;++i) co2.unsafe_indices()[i]=si[i]; for(int i=0;i<3;++i) co2.unsafe_offsets()[i]=sofs[i]; auto cso=cscout.edit().sparse_bytes();
    before=allocations; status=preprocess::log1p_into(sci,D::float64,L::csc,cso,D::float64,L::csc);
    if(!status || allocations!=before) return 18;
    if(!near(co2.values()[0],std::log(2)) || !near(co2.values()[1],std::log(3)) || !near(co2.values()[2],std::log(4)) || !near(co2.values()[3],std::log(7))) return 19;

    before=allocations; status=preprocess::normalize_into(in,D::float64,dst,D::float32,{true,0,false,.05},norm);
    if(status.code!=preprocess::QcCode::invalid_argument || allocations!=before) return 20;
    preprocess::NormalizationBuffers short_norm{factors,2,included,2,work,3,0};
    status=preprocess::normalization_factors(in,D::float64,automatic,short_norm);
    if(status.code!=preprocess::QcCode::insufficient_capacity || status.required!=3) return 21;
    status=preprocess::scale_into(cs,D::float64,L::csr,csout,D::float32,L::csr,{true,false,0,nullptr},mean,deviation,2,stored,2);
    if(status.code!=preprocess::QcCode::invalid_argument) return 22;

    double overlap_data[]={1,2,3,4,5,6}; auto overlap_in=anndatapp::Array::borrow_dense(overlap_data,{2,2},{2,1});
    double sentinels_zero[]={91,92,93,94}; auto output_zero=anndatapp::Array::borrow_dense(sentinels_zero,{2,2},{0,1});
    double sentinels_one[]={91,92,93,94}; auto output_one=anndatapp::Array::borrow_dense(sentinels_one,{2,2},{1,1});
    auto overlap_src=overlap_in.view().bytes(); auto oz=output_zero.edit().bytes(); auto oo=output_one.edit().bytes();
    before=allocations; status=preprocess::log1p_into(overlap_src,D::float64,oz,D::float64);
    if(status.code!=preprocess::QcCode::invalid_argument || allocations!=before || sentinels_zero[0]!=91 || sentinels_zero[1]!=92 || sentinels_zero[2]!=93 || sentinels_zero[3]!=94) return 40;
    status=preprocess::log1p_into(overlap_src,D::float64,oo,D::float64);
    if(status.code!=preprocess::QcCode::invalid_argument || allocations!=before || sentinels_one[0]!=91 || sentinels_one[1]!=92 || sentinels_one[2]!=93 || sentinels_one[3]!=94) return 41;
    double partial[]={1,2,3,4,5}; auto partial_in=anndatapp::Array::borrow_dense(partial,{2,2},{2,1});
    auto partial_out=anndatapp::Array::borrow_dense(partial+1,{2,2},{2,1}); auto partial_dst=partial_out.edit().bytes(); auto partial_src=partial_in.view().bytes();
    before=allocations;
    status=preprocess::log1p_into(partial_src,D::float64,partial_dst,D::float64);
    if(status.code!=preprocess::QcCode::invalid_argument || allocations!=before || partial[0]!=1 || partial[1]!=2 || partial[2]!=3 || partial[3]!=4 || partial[4]!=5) return 50;
    auto exact=anndatapp::Array::borrow_dense(overlap_data,{2,2},{2,1}); auto exact_view=exact.view().bytes(); auto exact_dst=exact.edit().bytes();
    before=allocations;
    status=preprocess::log1p_into(exact_view,D::float64,exact_dst,D::float64);
    if(!status || allocations!=before || !near(overlap_data[0],std::log1p(1)) || !near(overlap_data[3],std::log1p(4))) return 51;

    auto square_csr=anndatapp::Array::csr<double,int32_t,int32_t>(2,2,2);
    auto sqe=square_csr.edit().csr<double,int32_t,int32_t>();
    sqe.values()[0]=1; sqe.values()[1]=2; sqe.unsafe_indices()[0]=0; sqe.unsafe_indices()[1]=1;
    sqe.unsafe_offsets()[0]=0; sqe.unsafe_offsets()[1]=1; sqe.unsafe_offsets()[2]=2;
    auto square_csc=anndatapp::Array::csc<double,int32_t,int32_t>(2,2,2);
    auto sqo=square_csc.edit().csc<double,int32_t,int32_t>();
    sqo.values()[0]=81; sqo.values()[1]=82; sqo.unsafe_indices()[0]=0; sqo.unsafe_indices()[1]=1;
    sqo.unsafe_offsets()[0]=0; sqo.unsafe_offsets()[1]=1; sqo.unsafe_offsets()[2]=2;
    auto square_src=square_csr.view().sparse_bytes(); auto square_dst=square_csc.edit().sparse_bytes();
    status=preprocess::log1p_into(square_src,D::float64,L::csr,square_dst,D::float64,L::csc);
    if(status.code!=preprocess::QcCode::invalid_argument || sqo.values()[0]!=81 || sqo.values()[1]!=82) return 52;
    auto changed_csr=anndatapp::Array::csr<double,int32_t,int32_t>(2,2,2);
    auto changed=changed_csr.edit().csr<double,int32_t,int32_t>();
    changed.values()[0]=71; changed.values()[1]=72; changed.unsafe_indices()[0]=1; changed.unsafe_indices()[1]=0;
    changed.unsafe_offsets()[0]=0; changed.unsafe_offsets()[1]=1; changed.unsafe_offsets()[2]=2;
    auto changed_dst=changed_csr.edit().sparse_bytes();
    status=preprocess::log1p_into(square_src,D::float64,L::csr,changed_dst,D::float64,L::csr);
    if(status.code!=preprocess::QcCode::invalid_structure || changed.values()[0]!=71 || changed.values()[1]!=72) return 53;
    return 0;
}
