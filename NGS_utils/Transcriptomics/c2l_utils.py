def read_and_qc(sample_name, path=None):
    r""" This function reads the data for one 10X spatial experiment into the anndata object.
    It also calculates QC metrics. Modify this function if required by your workflow.
    
    :param sample_name: Name of the sample
    :param path: path to data
    """
    
    adata = sc.read_visium(path,
                           count_file='filtered_feature_bc_matrix.h5', load_images=True)
    adata.obs['sample'] = sample_name
    adata.var['SYMBOL'] = adata.var_names
    adata.var.rename(columns={'gene_ids': 'ENSEMBL'}, inplace=True)
    adata.var_names = adata.var['ENSEMBL']
    #adata.var.drop(columns='ENSEMBL', inplace=True)
    adata.obs["condition"] = infer_condition(sample_name)
    
    # Calculate QC metrics
    from scipy.sparse import csr_matrix
    adata.X = adata.X.toarray()
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    adata.X = csr_matrix(adata.X)
    adata.var['mt'] = [gene.startswith('mt-') for gene in adata.var['SYMBOL']]
    adata.obs['mt_frac'] = adata[:, adata.var['mt'].tolist()].X.sum(1).A.squeeze()/adata.obs['total_counts']

    # add sample name to obs names
    adata.obs_names = adata.obs["sample"] \
                          + '_' + adata.obs_names
    adata.obs.index.name = 'spot_id'
    return adata

def plot_visium_obs_grid(
    adata_vis,
    samples,                 # list of Visium library_ids you want as columns (e.g. ["wt1-10405", ...])
    obs_cols,                # list of adata_vis.obs column names to plot as rows
    sample_mapping,          # {library_id -> sample_code_in_obs}
    *,
    batch_key="sample",      # column in .obs that holds sample codes
    img_key="hires",
    cmap="magma",
    size=1.3,
    alpha_img=0.55,
    # scaling
    scale=None,  # None | "per_col_quantile" | "global_quantile" | "fixed"
    q_low=0.0,
    q_high=99.5,
    vmin_fixed=None,
    vmax_fixed=None,
    title_fn=None            # optional: function(ctcol) -> y-axis label
):
    """
    Generic grid plotter for any numeric columns in adata_vis.obs.
    - One column per Visium library in `samples`.
    - One row per entry in `obs_cols`.
    - Works for cell2location abundances or score_genes outputs alike.

    Scaling modes:
      - per_col_quantile: vmin/vmax computed per plotted column (robust to outliers).
      - global_quantile: vmin/vmax computed across *all samples* for each column.
      - fixed: use vmin_fixed/vmax_fixed for all plots
      - row_zero_one: normalize across the selected libraries only
        (visual comparison of where a cell type (or score) is enriched across slides.)
    """
    
    cols_for_plot = list(obs_cols)
    ylabels = list(obs_cols)
    
    # --- compute vmin/vmax per column ---
    vmin_map, vmax_map = {}, {}
    
    if scale == "fixed":
        if vmin_fixed is None or vmax_fixed is None:
            raise ValueError("scale='fixed' requires vmin_fixed and vmax_fixed.")
        for col in obs_cols:
            vmin_map[col] = float(vmin_fixed)
            vmax_map[col] = float(vmax_fixed)
    
    elif scale in {"global_quantile", "per_col_quantile"}:
        for col in obs_cols:
            vals = adata_vis.obs[col].astype(float).values
            vmin_map[col] = float(np.nanpercentile(vals, q_low))
            vmax_map[col] = float(np.nanpercentile(vals, q_high))
    
    elif scale == "row_zero_one":
        # compute row-wise min/max using *selected libraries only*
        sample_codes = {str(sample_mapping[s]) for s in samples}
        mask = adata_vis.obs[batch_key].astype(str).isin(sample_codes)
    
        new_cols = []
        for col in obs_cols:
            if col not in adata_vis.obs.columns:
                new_cols.append(col)  # will be caught later in plotting loop
                continue
            vals = adata_vis.obs.loc[mask, col].astype(float).values
            lo = float(np.nanpercentile(vals, q_low))
            hi = float(np.nanpercentile(vals, q_high))
            if not np.isfinite(hi) or hi <= lo:
                hi = lo + 1e-9
    
            # avoid double suffixing if already normalized
            norm_name = col if col.endswith("_n01") else f"{col}_n01"
            x = adata_vis.obs[col].astype(float).values
            x01 = np.clip((x - lo) / (hi - lo), 0, 1)
            adata_vis.obs[norm_name] = x01
    
            vmin_map[norm_name] = 0.0
            vmax_map[norm_name] = 1.0
            new_cols.append(norm_name)
    
        cols_for_plot = new_cols  # <-- use normalized names for plotting
    
    else:
        assert scale == None, "No scaling used"

    # --- prepare slices per library (keep only its own spatial block) ---
    slides = {}
    for lib_id in samples:
        if lib_id not in sample_mapping:
            raise KeyError(f"{lib_id!r} missing in sample_mapping.")
        sid = str(sample_mapping[lib_id])

        sub = adata_vis[adata_vis.obs[batch_key].astype(str) == sid].copy()
        
        if sub.n_obs == 0:
            warnings.warn(f"No spots for sample_code={sid} (library {lib_id}); will plot empty.")
        # ensure correct spatial entry
        if "spatial" not in adata_vis.uns or lib_id not in adata_vis.uns["spatial"]:
            # try fallback if user’s keys differ
            keys = list(adata_vis.uns.get("spatial", {}).keys())
            fallback = lib_id if lib_id in keys else (sid if sid in keys else (keys[0] if keys else None))
            if fallback is None:
                raise KeyError("No spatial metadata found in adata_vis.uns['spatial'].")
            lib_for_plot = fallback
        else:
            lib_for_plot = lib_id
        sub.uns["spatial"] = {lib_for_plot: adata_vis.uns["spatial"][lib_for_plot]}
        slides[lib_id] = (sub, lib_for_plot)

    # plot
    R, C = len(cols_for_plot), len(samples)
    fig, axes = plt.subplots(R, C, figsize=(4*C, 4.5*R), squeeze=False)
    
    for i, col_to_plot in enumerate(cols_for_plot):
        for j, lib_id in enumerate(samples):
            sub, lib_for_plot = slides[lib_id]
            if col_to_plot not in sub.obs.columns:
                axes[i, j].axis("off")
                axes[i, j].set_title(f"{lib_id}\n(missing '{col_to_plot}')", fontsize=9)
                continue
    
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=FutureWarning)
                if scale is None:
                    sc.pl.spatial(
                        sub,
                        color=col_to_plot,
                        cmap=cmap,
                        img_key=img_key,
                        library_id=lib_for_plot,
                        size=size,
                        alpha_img=alpha_img,
                        show=False,
                        ax=axes[i, j],
                    )
                else:
                    sc.pl.spatial(
                        sub,
                        color=col_to_plot,
                        cmap=cmap,
                        img_key=img_key,
                        library_id=lib_for_plot,
                        size=size,
                        alpha_img=alpha_img,
                        vmin=vmin_map.get(col_to_plot, None),
                        vmax=vmax_map.get(col_to_plot, None),
                        show=False,
                        ax=axes[i, j],
                    )
            if i == 0:
                axes[i, j].set_title(lib_id, fontsize=10)
    
        # left-side row label should show the original name
        axes[i, 0].set_ylabel(ylabels[i], fontsize=10)

    plt.tight_layout()
    return fig, axes

def compare_symbol_vs_ensembl(adata_ref, adata_vis, sym2ensemble, symbol_col="SYMBOL", ensembl_col="ENSEMBL"):
    """
    Compare overlaps between reference (single-cell) and visium datasets 
    using gene SYMBOLs vs Ensembl IDs.

    Parameters
    ----------
    adata_ref : AnnData
        Reference single-cell AnnData object (var_names assumed to be symbols).
    adata_vis : AnnData
        Visium AnnData object with SYMBOL and ENSEMBL columns in .var.
    sym2ensemble : dict-like
        Mapping from SYMBOL to Ensembl ID.
    symbol_col : str, default "SYMBOL"
        Column in adata_vis.var that stores gene symbols.
    ensembl_col : str, default "ENSEMBL"
        Column in adata_vis.var that stores Ensembl IDs.

    Returns
    -------
    result : dict
        Dictionary with overlap stats and collision info.
    """

    def strip_version(x): 
        return x.split(".", 1)[0] if isinstance(x, str) else x

    # symbol overlap
    sym_vis  = pd.Index(adata_vis.var[symbol_col].astype(str))
    sym_ref  = pd.Index(adata_ref.var_names.astype(str))
    sym_common = pd.Index(np.intersect1d(sym_vis, sym_ref))

    # map common symbols → Ensembl
    ens_from_sym_common = (
        pd.Index(sym_common.map(sym2ensemble))
        .dropna()
        .map(str)
        .unique()
    )

    # direct Ensembl overlap
    ens_vis = pd.Index(adata_vis.var[ensembl_col].map(strip_version)).dropna().map(str)
    ens_sc  = pd.Index(pd.Series(adata_ref.var_names, dtype=str).map(sym2ensemble)).dropna().map(str)
    ens_common_direct = pd.Index(np.intersect1d(ens_vis, ens_sc))

    # compare sets
    same = set(ens_from_sym_common) == set(ens_common_direct)

    # check collisions
    df_map = pd.DataFrame({"symbol": sym_common})
    df_map["ensembl"] = df_map["symbol"].map(sym2ensemble)
    collisions = (
        df_map.dropna()
              .groupby("ensembl")["symbol"]
              .nunique()
              .sort_values(ascending=False)
    )
    n_collide = (collisions > 1).sum()

    result = {
        "ensembl_from_symbols": set(ens_from_sym_common),
        "ensembl_direct": set(ens_common_direct),
        "same": same,
        "n_from_symbols": len(ens_from_sym_common),
        "n_direct": len(ens_common_direct),
        "extra_map": sorted(set(ens_from_sym_common) - set(ens_common_direct)),
        "extra_direct": sorted(set(ens_common_direct) - set(ens_from_sym_common)),
        "n_collisions": n_collide,
        "collisions": collisions[collisions > 1].head(10)
    }
    return result

def clear_c2l_reference_signatures(adata):
    """Remove RegressionModel export artifacts from a reference AnnData."""
    # varm signatures/statistics
    for k in [
        "means_per_cluster_mu_fg",
        "q05_per_cluster_mu_fg",
        "q95_per_cluster_mu_fg",
        "stds_per_cluster_mu_fg",
    ]:
        if k in adata.varm: adata.varm.pop(k, None)

    # sometimes extra keys with similar naming exist; remove any *_per_cluster_mu_fg
    for k in list(adata.varm.keys()):
        if k.endswith("_per_cluster_mu_fg"):
            adata.varm.pop(k, None)

    # model blob
    adata.uns.pop("mod", None)

    print("[clean] Reference signature artifacts removed.")

def clear_c2l_visium_posteriors(adata):
    for k in [
        "means_cell_abundance_w_sf",
        "q05_cell_abundance_w_sf",
        "q95_cell_abundance_w_sf",
        "cell_abundance_w_sf",
    ]:
        adata.obsm.pop(k, None)
    adata.uns.pop("mod", None)
    print("[clean] Cleared previous c2l posteriors from Visium AnnData.")
