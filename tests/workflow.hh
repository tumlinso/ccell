#pragma once

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
    auto stats = pp::feature_stats(s.view());
    if (stats.observations != 4 || stats.mean.view().dense<double>().data()[0] != 2.5 ||
        std::abs(stats.variance.view().dense<double>().data()[0] - 5.0/3) > 1e-12) return 1;
    v::PreparedInput input({s.view(), u.view(), w.view()});
    v::Workspace workspace(input);
    auto output = ad::Array::empty<double>({4, 1});
    v::FitOptions options;
    options.trimming = v::Trimming::none;
    auto fit = v::reference::steady_state_into(input, workspace, output.edit(), options);
    if (std::abs(fit.gamma.view().dense<double>().data()[0] - 2) > 1e-12) return 2;
    auto values = output.view().dense<double>();
    for (int i = 0; i < 4; ++i) if (std::abs(values.data()[i]) > 1e-12) return 3;
    return 0;
}
