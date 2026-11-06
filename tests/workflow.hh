#pragma once
#include <stddef.h>
#include <stdint.h>

int workflow() {
    namespace ad = ccell::anndata;
    namespace pp = ccell::preprocess;
    namespace v = ccell::velocity;
    auto s = ad::Array::empty<double>({4, 1});
    auto u = ad::Array::empty<double>({4, 1});
    auto w = ad::Array::zeros<double>({4, 4});
    auto* sp = s.edit().dense<double>().data();
    auto* up = u.edit().dense<double>().data();
    auto* wp = w.edit().dense<double>().data();
    for (int i = 0; i < 4; ++i) { sp[i] = i + 1; up[i] = 2*sp[i]; wp[5*i] = 1; }
    auto descriptor = s.view().bytes();
    double mean[1], variance[1];
    uint64_t stats_scratch[1];
    size_t observations = 0;
    auto status = pp::feature_stats(descriptor, ad::DType::float64, nullptr,
        mean, variance, 1, stats_scratch, 1, &observations);
    if (!status || observations != 4 || mean[0] != 2.5 ||
        __builtin_fabs(variance[0] - 5.0/3) > 1e-12) return 1;
    v::FitOptions options;
    options.trimming = v::Trimming::none;
    auto ud = u.view().bytes(), wd = w.view().bytes();
    v::Input input{{&descriptor, nullptr, ad::DType::float64},
                   {&ud, nullptr, ad::DType::float64},
                   {&wd, nullptr, ad::DType::float64}};
    v::PreparationSize preparation;
    if (!v::preparation_size(input, preparation)) return 2;
    size_t indices[128]; v::Entry entries[128];
    if (preparation.indices > 128 || preparation.entries > 128) return 3;
    v::PreparedInput prepared;
    if (!v::prepare(input, {indices, 128, entries, 128}, prepared)) return 4;
    v::WorkspaceSize requirement;
    if (!v::workspace_size<double>(prepared, {}, false, requirement)) return 5;
    double work[512]; uint8_t mask[16];
    if (requirement.scalars > 512 || requirement.masks > 16) return 6;
    v::Workspace<double> workspace;
    if (!v::bind_workspace(prepared, {}, false, work, 512, mask, 16, workspace)) return 7;
    double gamma[1], offset[1], r2[1], velocity[4];
    bool selected[1]; uint8_t fit_status[1];
    v::FitResult result{gamma, offset, r2, selected, fit_status, 1};
    if (!v::reference::steady_state_into(prepared, workspace, result, velocity, 4, 1, options)) return 8;
    if (__builtin_fabs(gamma[0] - 2) > 1e-12) return 9;
    for (int i = 0; i < 4; ++i)
        if (__builtin_fabs(velocity[i]) > 1e-12) return 10;
    return 0;
}
