from __future__ import annotations
import scanpy as sc, anndata as ad

from anndata.io import read_csv

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages

import os, sys
import os.path as path

import numpy as np
import pandas as pd

import warnings

from scipy import sparse
from scipy.io import mmread
from scipy.sparse import csc_matrix, csr_matrix

import cell2location, scvi, torch

from pathlib import Path
from typing import Optional, Tuple, Union, Dict

def get_directory_names(path_):
    return [d for d in os.listdir(path_) if os.path.isdir(os.path.join(path_, d))]

def read_barcode_norm_list_for_sample(csv_path: str, sample_label: str) -> set:
    df = pd.read_csv(csv_path)
    df['sample'] = df['sample'].astype(str)
    sel = df['sample'].str.upper() == str(sample_label).split("-")[0].upper()
    df = df.loc[sel]
    df["sample_id"] = str(sample_label)
    if df.empty:
        raise ValueError(f"No rows for sample='{sample_label}' in {csv_path}")
    col = 'barcode_norm' if 'barcode_norm' in df.columns else 'barcode'
    df["barcode_sample"] = df["sample_id"].astype(str).str.lower() + "_" + df[col]
    return set(df["barcode_sample"].astype(str).str.strip())


def infer_condition(sample_name: str) -> str:
    s = sample_name.upper()
    if s.startswith("WT"):  return "WT"
    if s.startswith("MGB"): return "MYOCD_mutant"
    return "unknown"

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
    scale=None,  # None | "row_zero_one" | "per_col_quantile" | "global_quantile" | "fixed"
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

def save_obs_grid_pdf(
    adata_vis, samples, obs_cols, sample_mapping, out_pdf, rows_per_page=3, **kwargs
):
    """Paginate long obs_cols lists into a multi-page PDF."""
    with PdfPages(out_pdf) as pdf:
        for i in range(0, len(obs_cols), rows_per_page):
            batch = obs_cols[i:i+rows_per_page]
            fig, _ = plot_visium_obs_grid(adata_vis, samples, batch, sample_mapping, **kwargs)
            fig.suptitle(f"{i+1}-{i+len(batch)} / {len(obs_cols)}", y=0.995, fontsize=10)
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
    print(f"Saved: {out_pdf}")


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


DFLike = Union[pd.DataFrame, Path, str]

def _load_q05_df(x: DFLike, strip_prefix: Optional[str] = "q05cell_abundance_w_sf_") -> pd.DataFrame:
    """Load a q05 abundance matrix from path or pass-through a DataFrame."""
    if isinstance(x, (str, Path)):
        x = Path(x)
        df = pd.read_csv(x, index_col=0, compression="gzip" if str(x).endswith(".gz") else None)
    else:
        df = x.copy()

    if strip_prefix:
        df.columns = [c.replace(strip_prefix, "") if isinstance(c, str) and c.startswith(strip_prefix) else c
                      for c in df.columns]
    df.index = df.index.astype(str)
    return df


def _exclude_mgb1(df: pd.DataFrame) -> pd.DataFrame:
    mask = ~df.index.str.split("_").str[0].str.lower().str.startswith("mgb1")
    return df.loc[mask].copy()


def run_abundance_enrichment(
    q05_whole_embryo: DFLike,
    q05_bladder_only: DFLike,
    outdir: Optional[Path] = None,
    *,
    strip_prefix: Optional[str] = "q05cell_abundance_w_sf_",
    exclude_mgb1_flag: bool = True,
    presence_mean_cutoff: float = 0.05,
    presence_nspots_cutoff: int = 3,
    presence_q05_cutoff: float = 0.25,
    enrichment_log2_cut: float = 1.0,      # >= 2x
    embryo_mean_floor: Optional[float] = None,  # e.g., 0.01 to avoid tiny-denominator inflation
    make_plots: bool = True,
    save_outputs: bool = True,
) -> Tuple[Dict[str, int | float | str], pd.DataFrame]:
    """
    Compute bladder enrichment from Round-1 q05 posteriors (whole-embryo vs bladder-only).

    - Enriched: log2(bladder/embryo) >= enrichment_log2_cut AND >= presence_nspots_cutoff bladder spots with q05 > presence_q05_cutoff
    - Present-but-shared: (mean_bladder > presence_mean_cutoff OR >= presence_nspots_cutoff spots above q05 threshold)
      but NOT enriched
    - Retain for next round: enriched OR present-but-shared

    Parameters
    ----------
    q05_whole_embryo, q05_bladder_only
        Path/str to CSV(.gz) or DataFrame with rows=spots, cols=cell states (q05).
    strip_prefix
        Remove this prefix from columns if present (default matches cell2location export).
    exclude_mgb1_flag
        Drop MGB1 spots by index prefix.
    embryo_mean_floor
        Optional floor on mean embryo abundance (e.g., 0.01) to avoid inflated ratios.
    make_plots, save_outputs
        Toggle plotting and saving. If `outdir` is None, saving is disabled.

    Returns
    -------
    summary : dict
    table   : pd.DataFrame (full results, sorted by log2_enrichment)
    """
    # I/O
    q05_emb = _load_q05_df(q05_whole_embryo, strip_prefix=strip_prefix)
    q05_bla = _load_q05_df(q05_bladder_only, strip_prefix=strip_prefix)

    if exclude_mgb1_flag:
        q05_emb = _exclude_mgb1(q05_emb)
        q05_bla = _exclude_mgb1(q05_bla)

    # Align columns
    common_cols = [c for c in q05_bla.columns if c in q05_emb.columns]
    if not common_cols:
        raise ValueError("No overlapping cell-state columns between embryo and bladder matrices.")
    q05_emb = q05_emb[common_cols].copy()
    q05_bla = q05_bla[common_cols].copy()

    # Means + enrichment
    eps = 1e-9
    mean_emb = q05_emb.mean(axis=0)
    mean_bla = q05_bla.mean(axis=0)

    if embryo_mean_floor is not None:
        mean_emb = mean_emb.clip(lower=float(embryo_mean_floor))

    enr = (mean_bla + eps) / (mean_emb + eps)
    log2_enr = np.log2(enr)

    # Presence support in bladder
    nspots_hi = (q05_bla > float(presence_q05_cutoff)).sum(axis=0)
    present_flag = (mean_bla > float(presence_mean_cutoff)) | (nspots_hi >= int(presence_nspots_cutoff))

    # Enrichment flag (with required bladder support)
    enriched_flag = (log2_enr >= float(enrichment_log2_cut)) & (nspots_hi >= int(presence_nspots_cutoff))

    # Retain either enriched OR clearly present (shared)
    retain_flag = enriched_flag | present_flag
    present_shared = present_flag & ~enriched_flag

    # Table
    table = (
        pd.DataFrame({
            "mean_q05_bladder": mean_bla,
            "mean_q05_embryo": mean_emb,
            "enrichment_ratio": enr,
            "log2_enrichment": log2_enr,
            f"bladder_spots_q05>{presence_q05_cutoff}": nspots_hi,
            "is_enriched": enriched_flag,
            "is_present_shared": present_shared,
            "retain_for_next_round": retain_flag,
        })
        .sort_values("log2_enrichment", ascending=False)
    )

    # Outputs
    if outdir is None:
        save_outputs = False
    if save_outputs:
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        table.to_csv(outdir / "bladder_enrichment_full.tsv", sep="\t")
        table[table["is_enriched"]].to_csv(outdir / "bladder_enrichment_only_enriched.tsv", sep="\t")

    # Plots
    if make_plots:
        # enriched-only barh
        enr_only = table[table["is_enriched"]]
        plt.figure(figsize=(9, max(3, 0.35 * len(enr_only))))
        xx = enr_only.sort_values("log2_enrichment", ascending=True)
        plt.barh(xx.index, xx["log2_enrichment"].values)
        plt.axvline(float(enrichment_log2_cut), linestyle="--")
        plt.xlabel("log2(bladder / embryo) q05")
        plt.title("Bladder-enriched lineages")
        plt.tight_layout()
        if save_outputs:
            plt.savefig(outdir / "bladder_enrichment_only_enriched.png", dpi=180)
            plt.close()

        # composition bar (enriched / shared / dropped)
        n_enr = int(table["is_enriched"].sum())
        n_shared = int(table["is_present_shared"].sum())
        n_drop = int((~table["retain_for_next_round"]).sum())
        plt.figure(figsize=(4, 4))
        plt.bar(["Enriched", "Shared", "Dropped"], [n_enr, n_shared, n_drop])
        plt.title("Lineage categories")
        plt.tight_layout()
        if save_outputs:
            plt.savefig(outdir / "lineage_categories.png", dpi=180)
            plt.close()

    summary = {
        "n_embryo_spots": int(q05_emb.shape[0]),
        "n_bladder_spots": int(q05_bla.shape[0]),
        "n_lineages": int(len(common_cols)),
        "n_enriched": int(table["is_enriched"].sum()),
        "n_shared_present": int(table["is_present_shared"].sum()),
        "n_dropped": int((~table["retain_for_next_round"]).sum()),
        "fig_enriched_bar": str(outdir / "bladder_enrichment_only_enriched.png") if save_outputs else "",
        "fig_category_bar": str(outdir / "lineage_categories.png") if save_outputs else "",
        "tsv_full": str(outdir / "bladder_enrichment_full.tsv") if save_outputs else "",
        "tsv_enriched_only": str(outdir / "bladder_enrichment_only_enriched.tsv") if save_outputs else "",
    }
    return summary, table

from cell2location.plt import plot_spatial

def make_slide_for_library(adata_vis, lib_id, sample_mapping, *, batch_key="sample", img_key="hires"):
    """
    Return a single-sample AnnData 'slide' configured so cell2location.plot_spatial uses the right image block.
    """
    if lib_id not in sample_mapping:
        raise KeyError(f"{lib_id!r} missing in sample_mapping.")
    sid = str(sample_mapping[lib_id])

    sub = adata_vis[adata_vis.obs[batch_key].astype(str) == sid].copy()
    if sub.n_obs == 0:
        warnings.warn(f"No spots for sample_code={sid} (library {lib_id}).")
    # ensure correct spatial entry
    if "spatial" not in adata_vis.uns or lib_id not in adata_vis.uns["spatial"]:
        # fallbacks if keys differ
        keys = list(adata_vis.uns.get("spatial", {}).keys())
        fallback = lib_id if lib_id in keys else (sid if sid in keys else (keys[0] if keys else None))
        if fallback is None:
            raise KeyError("No spatial metadata found in adata_vis.uns['spatial'].")
        lib_for_plot = fallback
    else:
        lib_for_plot = lib_id
    sub.uns["spatial"] = {lib_for_plot: adata_vis.uns["spatial"][lib_for_plot]}
    # pass through which image resolution to use
    if img_key is not None:
        sub.uns["spatial"][lib_for_plot]["use_quality"] = img_key
    return sub, lib_for_plot


def plot_visium_multicluster_per_sample(
    adata_vis,
    samples,
    clust_labels,                 # e.g. ['T_CD4+_naive','B_naive','FDC']  OR dict {lib_id: [..]}
    sample_mapping,               # {library_id -> sample_code_in_obs}
    *,
    batch_key="sample",
    img_key="hires",
    circle_diameter=6,
    max_color_quantile=0.992,
    colorbar_position="right",
    style="fast",
    show_img=True,
    return_figs=True
):
    """
    For each sample, overlay multiple cluster abundance columns on a single slide,
    using cell2location.plt.plot_spatial(color=[...], labels=[...]).

    clust_labels:
      - list[str] -> same set for all samples
      - dict[str, list[str]] -> per-sample customization
    Returns: dict {lib_id: (fig)} if return_figs else None
    """
    figs = {}
    # normalize clust_labels into per-sample dict
    if isinstance(clust_labels, dict):
        per_sample = {lib_id: [str(x) for x in clust_labels.get(lib_id, [])] for lib_id in samples}
    else:
        shared = [str(x) for x in clust_labels]
        per_sample = {lib_id: shared for lib_id in samples}

    for lib_id in samples:
        cols = per_sample[lib_id]
        if len(cols) == 0:
            warnings.warn(f"{lib_id}: no cluster columns provided; skipping.")
            continue

        # ensure columns exist in obs
        missing = [c for c in cols if c not in adata_vis.obs.columns]
        if missing:
            warnings.warn(f"{lib_id}: missing columns in .obs -> {missing}; they will be ignored.")
            cols = [c for c in cols if c in adata_vis.obs.columns]
            if not cols:
                continue

        slide, lib_for_plot = make_slide_for_library(
            adata_vis, lib_id, sample_mapping, batch_key=batch_key, img_key=img_key
        )

        # labels shown on the figure legend/colorbar (can differ from column keys)
        labels = [str(c) for c in cols]

        with plt.rc_context({"figure.figsize": (14, 14)}):
            fig = plot_spatial(
                adata=slide,
                color=cols,
                labels=labels,
                show_img=show_img,
                style=style,
                max_color_quantile=max_color_quantile,
                circle_diameter=circle_diameter,
                colorbar_position=colorbar_position,
            )
            fig.suptitle(f"{lib_id} — multi-cluster overlay", y=0.98, fontsize=12)
            if return_figs:
                figs[lib_id] = fig

    return figs if return_figs else None


from matplotlib.patches import Wedge

def plot_spatial_pies_c2lprop(
    adata,
    cols,                          # foreground cell types to color
    *,
    den_cols=None,                 # ALL c2l types for denominator (default: uns["mod"]["factor_names"])
    ax=None,
    fig=None,
    img_key="hires",
    show_img=True,
    img_alpha=0.9,
    include_other=True,
    other_label="Other",
    other_color=(0.7, 0.7, 0.7, 0.9),
    top_n=None,
    min_display_frac=0.0,
    rescale_shown=False,           # if True -> no Other
    radius_px=6,
    pad_px=40,
    crop=True,
    cmap_name="tab20",
    palette=None,                  # dict {celltype: color}
    legend=True,
    legend_ncols=2,
    style="fast",
    # --- labeling options (modular) ---
    label_from=None,               # None | "raw" | "drawn"
    label_min_frac=0.5,
    label_min_margin=0.10,
    label_min_total=None,
    write_labels=False,
    label_prefix="spot",           # writes e.g. spot_label, spot_label_frac, ...
):
    """
    Pies reflect true per-spot proportions:
        frac(ct) = abundance(ct) / sum(all_ct_abundances_in_that_spot)

    top_n + min_display_frac move mass into 'Other' (unless rescale_shown=True).
    Optionally compute labels using the same modular labeler from either:
        - 'raw'   : before any top_n/min_display_frac (preferred for analysis)
        - 'drawn' : after display trims (to mirror the figure)
    """
    if isinstance(cols, str):
        cols = [cols]
    cols = [c for c in cols if c in adata.obs.columns]
    if not cols:
        raise ValueError("None of the requested `cols` exist in adata.obs.")

    # denominator: all c2l types
    if den_cols is None:
        den_cols = list(adata.uns.get("mod", {}).get("factor_names", []))
        den_cols = [c for c in den_cols if c in adata.obs.columns]
        if not den_cols:
            den_cols = [c for c in adata.obs.columns if np.issubdtype(adata.obs[c].dtype, np.number)]

    # image + coords
    lib, block = next(iter(adata.uns["spatial"].items()))
    img = block.get("images", {}).get(img_key) if show_img else None
    sf = block["scalefactors"][f"tissue_{img_key}_scalef"]
    coords = adata.obsm["spatial"] * sf

    # --- build numerator and denominator (clip negatives) ---
    X_num = adata.obs[cols].to_numpy(float)
    X_num = np.nan_to_num(X_num, nan=0.0); X_num[X_num < 0] = 0.0

    X_den = adata.obs[den_cols].to_numpy(float)
    X_den = np.nan_to_num(X_den, nan=0.0); X_den[X_den < 0] = 0.0

    den = X_den.sum(axis=1, keepdims=True)
    den_safe = den.copy()
    den_safe[den_safe == 0] = 1.0

    # RAW proportions (science-faithful)
    F_raw = X_num / den_safe
    F_raw = np.clip(F_raw, 0, None)
    frac_raw_df = pd.DataFrame(F_raw, index=adata.obs_names, columns=cols)

    # --- display trims to get DRAWN fractions ---
    F = F_raw.copy()
    # Top-N within the requested cols
    if top_n is not None and 0 < top_n < F.shape[1]:
        kth = min(top_n, F.shape[1] - 1)
        keep_idx = np.argpartition(-F, kth=kth, axis=1)[:, :top_n]
        keep_mask = np.zeros_like(F, dtype=bool)
        keep_mask[np.arange(F.shape[0])[:, None], keep_idx] = True
    else:
        keep_mask = np.ones_like(F, dtype=bool)

    # Small slices to Other
    if min_display_frac and min_display_frac > 0:
        keep_mask &= (F >= float(min_display_frac))

    F_show = F * keep_mask
    shown_sum = F_show.sum(axis=1, keepdims=True)

    if rescale_shown:
        z = shown_sum.copy()
        z[z == 0] = 1.0
        F_draw = F_show / z
        other = np.zeros(F.shape[0], dtype=float)
    else:
        F_draw = F_show
        other = np.clip(1.0 - shown_sum.squeeze(), 0.0, 1.0)

    frac_drawn_df = pd.DataFrame(F_draw, index=adata.obs_names, columns=cols)

    # --- optional labeling (modular) ---
    labels_df = None
    if label_from in {"raw", "drawn"}:
        source_df = frac_raw_df if label_from == "raw" else frac_drawn_df
        totals_series = pd.Series(den.squeeze(), index=adata.obs_names) if label_min_total is not None else None
        labels_df = assign_labels_from_fractions(
            source_df,
            min_frac=label_min_frac,
            min_margin=label_min_margin,
            totals=totals_series,
            min_total=label_min_total,
            unassigned="Mixed/Low-confidence",
            return_details=True,
        )
        if write_labels:
            # write to obs with a clear prefix
            base = f"{label_prefix}_label"
            adata.obs[base] = labels_df["label"]
            adata.obs[f"{base}_frac"] = labels_df["label_frac"]
            adata.obs[f"{base}_second"] = labels_df["label_second"]
            adata.obs[f"{base}_margin"] = labels_df["label_margin"]

    # --- colors ---
    if palette is None:
        cmap = plt.colormaps.get_cmap(cmap_name)
        palette = {ct: cmap(i % cmap.N) for i, ct in enumerate(cols)}
    colors = [palette[c] for c in cols]

    # --- plotting ---
    with plt.style.context(style):
        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 10))
        elif fig is None:
            fig = ax.figure

        if img is not None:
            ax.imshow(img, origin="lower", alpha=img_alpha)

        if crop:
            xmin, ymin = coords.min(axis=0); xmax, ymax = coords.max(axis=0)
            ax.set_xlim(max(0, xmin - pad_px), xmax + pad_px)
            ax.set_ylim(max(0, ymin - pad_px), ymax + pad_px)
        ax.invert_yaxis()
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values(): sp.set_visible(False)

        for (x, y), row, oth in zip(coords, F_draw, other):
            start = 0.0
            # draw shown slices
            for j, frac in enumerate(row):
                f = float(frac)
                if f <= 0: 
                    continue
                theta = 360.0 * f
                ax.add_patch(Wedge((x, y), radius_px, start, start + theta,
                                   facecolor=colors[j], edgecolor="none", linewidth=0, alpha=0.95))
                start += theta
            # draw remainder as Other
            if include_other and (not rescale_shown) and (oth > 0):
                theta = 360.0 * float(oth)
                ax.add_patch(Wedge((x, y), radius_px, start, start + theta,
                                   facecolor=other_color, edgecolor="none", linewidth=0, alpha=0.95))

        if legend:
            from matplotlib.patches import Patch
            handles = [Patch(facecolor=palette[c], edgecolor="none", label=c) for c in cols]
            if include_other and not rescale_shown:
                handles.append(Patch(facecolor=other_color, edgecolor="none", label=other_label))
            ax.legend(handles=handles, title="Cell states", ncols=legend_ncols,
                      bbox_to_anchor=(1.02, 1.0), loc="upper left", frameon=False)

        fig.tight_layout()

    return fig, ax, frac_drawn_df, labels_df
