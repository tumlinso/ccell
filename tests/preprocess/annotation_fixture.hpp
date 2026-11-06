#pragma once

// Development-only adapter: the reference routines borrow buffers, while the
// qualification protocol exercises explicit AnnData ownership and annotations.
inline void run_annotations(anndatapp::Array& a, int mode, std::size_t n, std::size_t g,
                                       const std::vector<int>& groups) {
    auto data = anndatapp::AnnData::empty(n, g);
    data.set_x(std::move(a));
    std::vector<std::string> obs_names_order(n), var_names_order(g);
    for (std::size_t i = 0; i < n; ++i) obs_names_order[i] = "sample_" + std::to_string(n-i);
    for (std::size_t j = 0; j < g; ++j) var_names_order[j] = "gene_" + std::to_string(g-j);
    data.set_obs(anndatapp::Table::with_index(obs_names_order));
    data.set_var(anndatapp::Table::with_index(var_names_order));
    auto attach = [&](bool obs, const std::string& name, const auto& values) {
        using V = typename std::decay_t<decltype(values)>::value_type;
        auto arr = anndatapp::Array::empty<V>({values.size()});
        auto* dst = arr.edit().template dense<V>().data();
        for (std::size_t i = 0; i < values.size(); ++i) dst[i] = values[i];
        if (obs) data.set_obs_column(name, anndatapp::Column::numeric(std::move(arr)));
        else data.set_var_column(name, anndatapp::Column::numeric(std::move(arr)));
    };
    std::vector<std::int64_t> obs_sentinel(n), var_sentinel(g);
    for (std::size_t i = 0; i < n; ++i) obs_sentinel[i] = std::int64_t(i);
    for (std::size_t j = 0; j < g; ++j) var_sentinel[j] = std::int64_t(j);
    attach(true, "sentinel", obs_sentinel); attach(false, "sentinel", var_sentinel);
    data.set_layer("counts", data.x().copy());

    std::vector<double> ct(n), gt(g), mean(g), dropout(g), lct(n), lcd(n), lgt(g), lgm(g),
        group_total(n), group_pct(n), group_log(n), top(n*3), scratch(g);
    std::vector<std::uint64_t> cd(n), gd(g);
    std::vector<std::uint8_t> mask(g);
    for (std::size_t j = 0; j < g; ++j) mask[j] = groups[j] != 0;
    const std::size_t tops[] = {2, 4, 8};
    preprocess::QcOptions qo{mask.data(), 1, tops, 3};
    preprocess::QcBuffers qb;
    qb.cell_total = ct.data(); qb.cell_detected = cd.data(); qb.cells_capacity = n;
    qb.gene_total = gt.data(); qb.gene_detected = gd.data(); qb.genes_capacity = g;
    qb.mean = mean.data(); qb.dropout = dropout.data(); qb.log_cell_total = lct.data(); qb.log_cell_detected = lcd.data();
    qb.log_gene_total = lgt.data(); qb.log_gene_mean = lgm.data();
    qb.group_total = group_total.data(); qb.group_percent = group_pct.data(); qb.group_log_total = group_log.data(); qb.group_capacity = n;
    qb.top_percent = top.data(); qb.top_capacity = top.size();
    if (!preprocess::qc(data.x().bytes(), data.x().dtype(), qo, qb, scratch.data(), scratch.size()))
        throw std::runtime_error("annotation QC failed");
    attach(true,"total_counts",ct); attach(true,"n_genes_by_counts",cd);
    attach(true,"log1p_total_counts",lct); attach(true,"log1p_n_genes_by_counts",lcd);
    attach(true,"total_counts_control",group_total); attach(true,"pct_counts_control",group_pct);
    attach(true,"log1p_total_counts_control",group_log);
    for (std::size_t k = 0; k < 3; ++k) {
        std::vector<double> column(n);
        for (std::size_t i = 0; i < n; ++i) column[i] = top[i*3+k];
        attach(true,"pct_counts_in_top_"+std::to_string(tops[k])+"_genes",column);
    }
    attach(false,"total_counts",gt); attach(false,"n_cells_by_counts",gd);
    attach(false,"mean_counts",mean); attach(false,"pct_dropout_by_counts",dropout);
    attach(false,"log1p_total_counts",lgt); attach(false,"log1p_mean_counts",lgm);

    std::vector<double> factors(n), norm_scratch(n);
    std::vector<std::uint8_t> included(g);
    preprocess::NormalizationOptions no; no.target_is_set = true; no.target_sum = 100;
    preprocess::NormalizationBuffers nb{factors.data(), n, included.data(), g, norm_scratch.data(), n};
    if (!preprocess::normalize_into(data.x().bytes(), data.x().dtype(), data.edit_x().bytes(), data.x().dtype(), no, nb))
        throw std::runtime_error("annotation normalization failed");
    attach(true,"size_factors",factors);
    const double base = mode ? (mode == 1 ? 2.0 : .5) : 0;
    if (!preprocess::log1p_into(data.x().bytes(), data.x().dtype(), data.edit_x().bytes(), data.x().dtype(), base))
        throw std::runtime_error("annotation log1p failed");
    auto log_meta = anndatapp::Element::mapping();
    log_meta.set("base", mode ? anndatapp::Element(base) : anndatapp::Element(nullptr));
    data.set_uns("log1p", std::move(log_meta));

    std::vector<double> hvmean(g), dispersion(g), normdisp(g);
    std::vector<std::uint8_t> selected(g);
    preprocess::HvgOptions ho; ho.has_top = true; ho.n_top_genes = 12; ho.n_bins = 10;
    ho.has_log_base = mode != 0; ho.log_base = base;
    preprocess::HvgOutput hout; hout.mean = hvmean.data(); hout.dispersion = dispersion.data();
    hout.normalized_dispersion = normdisp.data(); hout.selected = selected.data(); hout.capacity = g;
    std::size_t hbytes = 0;
    if (preprocess::hvg_scratch_bytes(g, 1, ho, hbytes).code != preprocess::HvgCode::ok)
        throw std::runtime_error("annotation HVG scratch failed");
    qualification::Scratch hscratch(hbytes);
    auto x = data.x().bytes();
    if (preprocess::hvg_dense(x.data(), data.x().dtype(), n, g, x.stride_bytes(0), x.stride_bytes(1), nullptr, 1,
                              ho, hout, hscratch.data, hbytes).code != preprocess::HvgCode::ok)
        throw std::runtime_error("annotation HVG failed");
    attach(false,"means",hvmean); attach(false,"dispersions",dispersion); attach(false,"dispersions_norm",normdisp);
    auto bools = anndatapp::Array::empty<bool>({g});
    auto* bool_data = bools.edit().dense<bool>().data();
    for (std::size_t j = 0; j < g; ++j) bool_data[j] = selected[j] != 0;
    data.set_var_column("highly_variable", anndatapp::Column::numeric(std::move(bools)));
    auto hvg_meta = anndatapp::Element::mapping(); hvg_meta.set("flavor", anndatapp::Element("seurat"));
    data.set_uns("hvg", std::move(hvg_meta));

    std::vector<double> scale_mean(g), stddev(g);
    std::vector<std::uint64_t> scale_scratch(g);
    if (!preprocess::scale_into(data.x().bytes(), data.x().dtype(), data.edit_x().bytes(), data.x().dtype(), {},
                                scale_mean.data(), stddev.data(), g, scale_scratch.data(), g))
        throw std::runtime_error("annotation scaling failed");
    attach(false,"mean",scale_mean); attach(false,"std",stddev);
    const char* obs_names[] = {"total_counts", "n_genes_by_counts", "log1p_total_counts",
        "log1p_n_genes_by_counts", "total_counts_control", "pct_counts_control",
        "log1p_total_counts_control", "pct_counts_in_top_2_genes", "pct_counts_in_top_4_genes",
        "pct_counts_in_top_8_genes", "size_factors"};
    const char* var_names[] = {"total_counts", "n_cells_by_counts", "mean_counts", "pct_dropout_by_counts",
        "log1p_total_counts", "log1p_mean_counts", "means", "dispersions", "dispersions_norm", "highly_variable", "mean", "std"};
    for (auto name : obs_names) { std::string key = std::string("obs_") + name; emit(key.c_str(), data.obs_column(name).values()); }
    for (auto name : var_names) { std::string key = std::string("var_") + name; emit(key.c_str(), data.var_column(name).values()); }
    auto uns = data.uns();
    bool log_ok = uns.contains("log1p") && uns.at("log1p").contains("base");
    if (log_ok) { auto b = uns.at("log1p").at("base");
        log_ok = mode ? b.kind() == anndatapp::ElementKind::floating && b.floating() == base : b.kind() == anndatapp::ElementKind::null; }
    const bool hvg_ok = uns.contains("hvg") && uns.at("hvg").contains("flavor") && uns.at("hvg").at("flavor").string() == "seurat";
    std::cout << "metadata 1 2 " << log_ok << ' ' << hvg_ok << '\n';
    bool obs_aligned = true, var_aligned = true;
    for (std::size_t i = 0; i < n; ++i) obs_aligned &= data.obs().name(i) == obs_names_order[i];
    for (std::size_t j = 0; j < g; ++j) var_aligned &= data.var().name(j) == var_names_order[j];
    std::cout << "axis_alignment 1 2 " << obs_aligned << ' ' << var_aligned << '\n';
    emit("obs_sentinel", data.obs_column("sentinel").values()); emit("var_sentinel", data.var_column("sentinel").values());
    emit("layer_counts", data.layer("counts")); emit_matrix(data.x());
}
