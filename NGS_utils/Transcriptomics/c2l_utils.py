from __future__ import annotations
import scanpy as sc, anndata as ad

from anndata.io import read_csv

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm, colors
from matplotlib.patches import Wedge
from matplotlib.colors import to_rgba

import re
from datetime import datetime

import os, sys
import os.path as path

import numpy as np
import pandas as pd

import warnings

from scipy import sparse
from scipy.io import mmread
from scipy.sparse import csc_matrix, csr_matrix

import cell2location, scvi, torch
from cell2location.plt import plot_spatial

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

# ---------- internals ----------
_PREFIX_PAT = re.compile(r'^(?:means|mean|q05|q95)_?cell_abundance_w_sf_')

def _strip_c2l_prefixes(cols: pd.Index) -> pd.Index:
    return pd.Index([_PREFIX_PAT.sub('', str(c)) for c in cols])

def _ensure_df_from_obsm(adata, obsm_key: str, factor_names: list) -> pd.DataFrame:
    X = adata.obsm[obsm_key]
    if isinstance(X, pd.DataFrame):
        X = X.reindex(index=adata.obs_names)
        X.columns = _strip_c2l_prefixes(pd.Index(X.columns))
        if set(factor_names).issubset(X.columns):
            X = X.loc[:, factor_names]
        elif X.shape[1] == len(factor_names):
            X = X.copy(); X.columns = factor_names
        else:
            missing = [f for f in factor_names if f not in X.columns]
            raise ValueError(f"{obsm_key}: columns don’t match factor_names; missing e.g. {missing[:10]}")
        return X.astype(float)
    else:
        return pd.DataFrame(X, index=adata.obs_names, columns=factor_names).astype(float)

def _available_c2l_keys(adata):
    keys = {}
    if "means_cell_abundance_w_sf" in adata.obsm: keys["mean"] = "means_cell_abundance_w_sf"
    if "cell_abundance_w_sf"      in adata.obsm and "mean" not in keys:
        keys["mean"] = "cell_abundance_w_sf"  # legacy
    if "q05_cell_abundance_w_sf"  in adata.obsm: keys["q05"] = "q05_cell_abundance_w_sf"
    if "q95_cell_abundance_w_sf"  in adata.obsm: keys["q95"] = "q95_cell_abundance_w_sf"
    return keys

# ---------- main (improved) ----------
def aggregate_c2l_celltypes(
    adata,
    map_df,
    *,
    source_col="celltype_sub",
    desired_col="c2l_celltype_sub_short",
    flavors=("q05","mean","q95"),          # which posterior flavors to aggregate if available
    prefix="c2l_",                         # prefix for output .obsm keys
    write_to_obs=False,                    # default OFF (safer)
    obs_flavor=None,                       # if writing to .obs, which single flavor to use
    obs_suffix=False,                      # add '__{flavor}' suffix when writing to .obs
    registry_key="c2l_group_registry"      # where to log provenance in .uns
):
    """
    Aggregate c2l posteriors across fine factors into grouped cell types.

    Writes to .obsm (one matrix per available flavor):
        '{prefix}{flavor}_cell_abundance_w_sf'  e.g., 'c2l_q05_cell_abundance_w_sf'
        columns = grouped names from `desired_col` (UNprefixed; readers can add 'c2l_' if desired)

    Optionally writes one chosen flavor to .obs (NOT recommended unless you really need it).

    Returns a dict with: groups mapping, written .obsm keys, and (optionally) written .obs columns.
    """
    # 1) inputs & mapping
    factor_names = list(adata.uns.get("mod", {}).get("factor_names", []))
    if not factor_names:
        raise ValueError("Missing adata.uns['mod']['factor_names'].")

    m = (map_df[[source_col, desired_col]]
         .dropna()
         .astype(str)
         .drop_duplicates(subset=[source_col]))
    groups = m.groupby(desired_col)[source_col].apply(list).to_dict()
    if not groups:
        raise ValueError("No groups found from map_df; check source_col / desired_col.")

    # 2) which posteriors exist
    avail = _available_c2l_keys(adata)
    if not avail:
        raise ValueError("No c2l matrices in .obsm (expected mean/cell_abundance_w_sf and/or q05/q95).")

    # limit to requested flavors that are present
    to_do = [(f, avail[f]) for f in flavors if f in avail]
    if not to_do:
        raise ValueError(f"None of requested flavors {flavors} exist. Available: {list(avail)}")

    obsm_written = []
    obs_cols_written = []

    # 3) aggregate for each flavor
    for flavor, src_key in to_do:
        X = _ensure_df_from_obsm(adata, src_key, factor_names)  # [spots × factors]

        agg = {}
        for new_name, old_list in groups.items():
            present = [c for c in old_list if c in X.columns]
            if present:
                agg[new_name] = X[present].sum(axis=1)
        if not agg:
            continue

        agg_df = pd.DataFrame(agg, index=adata.obs_names).astype(float)
        obsm_out = f"{prefix}{flavor}_cell_abundance_w_sf"   # e.g., 'c2l_q05_cell_abundance_w_sf'
        adata.obsm[obsm_out] = agg_df
        obsm_written.append(obsm_out)

    # 4) optional .obs write (single flavor only, discouraged by default)
    if write_to_obs:
        if obs_flavor is None:
            # pick the first actually written flavor
            obs_flavor = next((f for f, _ in to_do), None)
        if obs_flavor is None:
            raise ValueError("Requested write_to_obs but no flavor available to write.")
        key = f"{prefix}{obs_flavor}_cell_abundance_w_sf"
        if key not in adata.obsm:
            raise ValueError(f"Requested obs_flavor='{obs_flavor}' not aggregated into .obsm.")
        df = adata.obsm[key]
        cols = df.columns
        out_cols = (cols if not obs_suffix
                    else pd.Index([f"{c}__{obs_flavor}" for c in cols]))
        # write with prefix in obs so they look like 'c2l_<group>'
        obs_block = df.copy()
        obs_block.columns = pd.Index([f"{prefix}{c}" for c in out_cols])
        # atomic replace-insert
        tmp = adata.obs.copy()
        tmp[obs_block.columns] = obs_block.reindex(tmp.index)
        adata.obs = tmp
        obs_cols_written = list(obs_block.columns)

    # 5) provenance
    log = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "source_keys": {f: avail[f] for f, _ in to_do},
        "obsm_written": obsm_written,
        "wrote_obs": write_to_obs,
        "obs_flavor": obs_flavor if write_to_obs else None,
        "mapping_cols": {"source_col": source_col, "desired_col": desired_col},
        "n_groups": len(groups),
    }
    adata.uns.setdefault(registry_key, []).append(log)

    return {
        "groups": groups,
        "obsm_keys_written": obsm_written,
        "obs_flavor_written": obs_flavor if write_to_obs else None,
        "obs_cols_written": obs_cols_written,
        "log": log,
    }

def build_grouped_palette(cols_to_show, base_colors, fallback_cmap="tab20"):
    """
    cols_to_show: e.g. ['c2l_Urothelial Cells','c2l_Fibroblasts1', ...]
    base_colors:  dict keyed by group names WITHOUT 'c2l_' (spaces OK)
    Returns dict keyed exactly by cols_to_show.
    """
    # normalize keys in base map (case/underscore/space agnostic)
    def norm(s): return re.sub(r'\s+', '_', s.strip().lower())
    norm_base = {norm(k): v for k, v in base_colors.items()}

    cmap = plt.colormaps.get_cmap(fallback_cmap)
    pal = {}
    auto_i = 0
    for col in cols_to_show:
        raw = col.removeprefix("c2l_")
        candidates = [raw, raw.replace('_',' '), raw.replace(' ','_')]
        hit = None
        for cand in candidates:
            k = norm(cand)
            if k in norm_base:
                hit = norm_base[k]; break
        if hit is None:
            pal[col] = cmap(auto_i % cmap.N); auto_i += 1  # auto color if not provided
        else:
            pal[col] = to_rgba(hit)
    return pal

def plot_spatial_pies_c2lprop(
    adata,
    cols,                          # foreground cell types to color
    *,
    den_cols=None,                 # ALL c2l types for denominator (default: uns["mod"]["factor_names"])
    obsm_key="c2l_q05_cell_abundance_w_sf",
    den_from="post",                 # "post" (default) or "raw"
    raw_obsm_key="q05_cell_abundance_w_sf",
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
    label_prefix="spot",          # writes e.g. spot_label, spot_label_frac, ...,
    savepath=None,
    savefile=None
):
    """
    Pies reflect true per-spot proportions:
        frac(ct) = abundance(ct) / sum(all_ct_abundances_in_that_spot)

    top_n + min_display_frac move mass into 'Other' (unless rescale_shown=True).
    Optionally compute labels using the same modular labeler from either:
        - 'raw'   : before any top_n/min_display_frac (preferred for analysis)
        - 'drawn' : after display trims (to mirror the figure)
    """
    
    POST = adata.obsm[obsm_key]
    if not isinstance(POST, pd.DataFrame):
        POST = pd.DataFrame(POST, index=adata.obs_names)

    # --- denominator source ---
    if den_from == "post":
        if den_cols is None:
            den_cols = list(POST.columns)
        else:
            missing = [c for c in den_cols if c not in POST.columns]
            if missing:
                raise KeyError(f"den_cols not in {obsm_key}: {missing}")
        X_den = POST[den_cols].to_numpy(float)

    elif den_from == "raw":
        RAW = adata.obsm[raw_obsm_key]
        if not isinstance(RAW, pd.DataFrame):
            RAW = pd.DataFrame(RAW, index=adata.obs_names)
        # strip standard c2l prefix once for raw matrices
        import re
        _pat = re.compile(r'^(?:means|mean|q05|q95)_?cell_abundance_w_sf_')
        RAW = RAW.copy()
        RAW.columns = pd.Index([_pat.sub('', str(c)) for c in RAW.columns])

        use_den = den_cols or list(RAW.columns)
        miss = [c for c in use_den if c not in RAW.columns]
        if miss:
            raise KeyError(f"den_cols not in {raw_obsm_key}: {miss}")
        X_den = RAW[use_den].to_numpy(float)
    else:
        raise ValueError("den_from must be 'post' or 'raw'")

    # ---- image + coordinates ----
    lib, block = next(iter(adata.uns["spatial"].items()))
    img = block.get("images", {}).get(img_key) if show_img else None
    sf = block["scalefactors"][f"tissue_{img_key}_scalef"]
    coords = adata.obsm["spatial"] * sf

    # --- build numerator and denominator (clip negatives) ---
    X_num = POST[cols].to_numpy(float)
    X_num = np.nan_to_num(X_num, nan=0.0); X_num[X_num < 0] = 0.0

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
        if savefile is not None and savepath is not None:
            plt.savefig(savepath / savefile, dpi=300)
            plt.close()
        
    return fig, ax, frac_drawn_df, labels_df

def c2l_abundance_by_compartment_genotype_sectionwise(
    adata,
    meta_data: pd.DataFrame,           # columns: sample_id, barcode_norm, region
    *,
    cols=None,                         # list of 'c2l_*' names to include; if None, auto from source
    flavor="q05",                      # 'q05' | 'mean' | 'q95'
    ensure_flavor=False,               # auto-build .obsm flavor if missing (no .obs writes)
    map_df=None,                       # required if ensure_flavor=True
    aggregate_fn=None,                 # function handle to aggregate_c2l_groups_refined
    sample_id_col="sample_id",
    barcode_col="barcode_norm",
    region_col="region",
    compartments_order=None,           # e.g. ["Outer","Middle","Inner","Umblical"]
    genotype_key="genotype",
    genotype_order=None,               # e.g. ["WT","mutant"]
    section_key="section_id",
    exclude_sections=None
):
    """
    Section-equalized absolute abundance per (compartment, genotype).

    Steps:
      1) Choose aggregated c2l matrix (prefer .obsm 'c2l_{flavor}_cell_abundance_w_sf').
      2) Map each spot to a compartment (region) via meta_data (sample_id + barcode_norm).
      3) For EACH section_id: mean abundance per (compartment, cell).
      4) For EACH genotype: average the section-level tables (equal weight per section).

    Returns:
      tables_by_genotype: dict[genotype] -> DataFrame (rows=cells, cols=compartments)
      per_section_tables: dict[section_id] -> DataFrame (rows=cells, cols=compartments)
      used_source: string describing which source matrix was used
    """

    # --- optionally ensure requested flavor exists in .obsm ---
    obsm_key_pref = f"c2l_{flavor}_cell_abundance_w_sf"
    if ensure_flavor and obsm_key_pref not in adata.obsm:
        if aggregate_fn is None or map_df is None:
            raise ValueError("ensure_flavor=True requires aggregate_fn and map_df.")
        aggregate_fn(
            adata, map_df, prefix="c2l_",
            write_to_obs=False,  # do not touch .obs
            obs_from=None
        )

    # --- choose source (prefer requested .obsm, else any c2l_* in .obsm, else .obs c2l_*) ---
    if obsm_key_pref not in adata.obsm:
        raise KeyError(f"Required matrix '{obsm_key_pref}' not found in .obsm. "
                       f"Set ensure_flavor=True or run the aggregator to create it.")
    POST, used_source = adata.obsm[obsm_key_pref], obsm_key_pref

    if not isinstance(POST, pd.DataFrame):
        POST = pd.DataFrame(POST, index=adata.obs_names)

    # subset/order requested cell IDs
    if cols is None:
        use_cols = list(POST.columns)
    else:
        use_cols = [c for c in cols if c in POST.columns]
        if not use_cols:
            raise ValueError(f"None of requested `cols` found in source ({used_source}).")
    POST = POST[use_cols].clip(lower=0)

    # optionally drop excluded sections
    obs_sub = adata.obs
    if exclude_sections:
        keep = ~obs_sub[section_key].astype(str).isin(exclude_sections)
        POST = POST.loc[keep]
        obs_sub = obs_sub.loc[keep]

    # ensure required obs keys
    for k in (section_key, genotype_key):
        if k not in obs_sub.columns:
            raise KeyError(f"adata.obs missing '{k}'")

    # --- map spots → region via composite key (sample_id + sep + barcode_norm) ---
    def _best_sep(F_idx: pd.Index, msub: pd.DataFrame):
        seps = ("_","-","","|",":")
        fset = set(F_idx.astype(str))
        best, hits = None, -1
        for s in seps:
            comp = (msub[sample_id_col].astype(str) + s + msub[barcode_col].astype(str)).values
            h = np.sum(np.isin(comp, list(fset)))
            if h > hits:
                best, hits = s, h
        if hits == 0:
            raise ValueError("Meta rows did not match spot indices; check sample_id/barcode_norm.")
        return best

    present_sections = set(obs_sub[section_key].astype(str))
    m = meta_data[meta_data[sample_id_col].astype(str).isin(present_sections)].copy()
    sep = _best_sep(POST.index, m)
    m["__key"] = m[sample_id_col].astype(str) + sep + m[barcode_col].astype(str)
    lut = m.set_index("__key")[region_col].astype(str)

    regions = pd.Series(POST.index.astype(str), index=POST.index).map(lut)
    mask = regions.notna()
    if not mask.any():
        raise ValueError("Regions could not be mapped. Check meta_data join columns.")
    POST = POST.loc[mask]
    regions = regions.loc[mask]

    # pull section + genotype aligned to POST rows
    sections = obs_sub.loc[POST.index, section_key].astype(str)
    genos    = obs_sub.loc[POST.index, genotype_key].astype(str)

    # determine column (compartment) order
    if compartments_order is None:
        compartments_order = list(pd.unique(regions))

    # --- PER-SECTION means: for each section, average spots inside each region ---
    per_section_tables = {}  # section_id -> DataFrame (rows=cells, cols=compartments)
    for sec, idx in sections.groupby(sections).groups.items():
        F_sec = POST.loc[idx]
        reg_sec = regions.loc[idx]

        # group by region, mean across spots, then transpose to cells×regions
        reg_comp = (F_sec.assign(_reg=reg_sec.values)
                          .groupby("_reg").mean(numeric_only=True)
                          .T)  # rows=cells, cols=regions

        # enforce dense, ordered columns
        for c in compartments_order:
            if c not in reg_comp.columns:
                reg_comp[c] = 0.0
        reg_comp = reg_comp[compartments_order]
        # enforce full row set/order
        reg_comp = reg_comp.reindex(use_cols).fillna(0.0)

        per_section_tables[sec] = reg_comp

    # --- EQUAL-WEIGHT across sections WITHIN genotype ---
    if genotype_order is None:
        genotype_order = list(pd.unique(genos))

    # collect list of sections per genotype
    sec2geno = sections.groupby(sections).agg(lambda s: genos.loc[s.index[0]])
    geno_to_secs = {g: [s for s, gg in sec2geno.items() if gg == g] for g in genotype_order}

    tables_by_genotype = {}
    for g in genotype_order:
        secs = geno_to_secs.get(g, [])
        if not secs:
            continue
        # simple (unweighted) mean of the section tables
        stack = np.stack([per_section_tables[s].to_numpy(float) for s in secs], axis=2)  # cells × comps × n_sec
        mean_gc = np.nanmean(stack, axis=2)
        df_g = pd.DataFrame(mean_gc, index=use_cols, columns=compartments_order)
        tables_by_genotype[g] = df_g

    return tables_by_genotype, per_section_tables, used_source


def plot_genotype_compartment_heatmaps(
    tables_by_genotype: dict,
    *,
    strip_prefix="c2l_",
    cmap="viridis",
    zscore_rows=False,          # False recommended for absolute comparability
    figsize_per=(6, 5),
    suptitle="Cell2location abundance per compartment (section-equalized)"
):
    # determine shared vmin/vmax across genotypes for fair comparison
    mats = [df.to_numpy(float) for df in tables_by_genotype.values()]
    big = np.hstack(mats) if mats else np.zeros((1,1))
    if zscore_rows:
        vmin = vmax = None
    else:
        vmax = np.nanmax(big) if np.isfinite(big).any() else 1.0
        vmin = 0.0

    n = len(tables_by_genotype)
    fig, axes = plt.subplots(1, n, figsize=(figsize_per[0]*n, figsize_per[1]), squeeze=False)
    axes = axes[0]

    for ax, (geno, df) in zip(axes, tables_by_genotype.items()):
        disp = df.copy()
        if strip_prefix:
            disp.index = [i.replace(strip_prefix, "") for i in disp.index]
        data = disp.to_numpy(float)

        if zscore_rows:
            mu = np.nanmean(data, axis=1, keepdims=True)
            sd = np.nanstd(data, axis=1, ddof=1, keepdims=True); sd[sd==0]=1.0
            data = (data - mu) / sd
            vmax_local = np.nanmax(np.abs(data)) if np.isfinite(data).any() else 1.0
            norm = colors.TwoSlopeNorm(vmin=-vmax_local, vcenter=0.0, vmax=vmax_local)
            cmap_use = "bwr"
        else:
            norm = None
            cmap_use = cmap

        im = ax.imshow(data, aspect="auto", cmap=cmap_use, norm=norm, vmin=None if zscore_rows else vmin, vmax=None if zscore_rows else vmax)
        ax.set_title(geno)
        ax.set_xlabel("Compartment")
        ax.set_xticks(np.arange(disp.shape[1])); ax.set_xticklabels(disp.columns, rotation=0)
        ax.set_yticks(np.arange(disp.shape[0])); ax.set_yticklabels(disp.index)


    # shared colorbar
    cax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    if zscore_rows:
        sm = axes[-1].images[-1]
        fig.colorbar(sm, cax=cax, label="Row z-score")
    else:
        import matplotlib.cm as cm
        sm = cm.ScalarMappable(norm=colors.Normalize(vmin=vmin, vmax=vmax), cmap=cmap)
        sm.set_array([])
        fig.colorbar(sm, cax=cax, label="Posterior abundance (mean)")

    fig.suptitle(suptitle, y=0.98)
    fig.tight_layout(rect=[0, 0, 0.9, 0.96])
    return fig, axes

# Δ (mutant − WT) from section-equalized tables
def plot_delta_heatmap_from_genotype_tables(tables_by_geno, *,
                                            wt_label="WT", mut_label="mutant",
                                            strip_prefix="c2l_", cmap="bwr",
                                            figsize=(8,6), title=None):
    delta = tables_by_geno[mut_label] - tables_by_geno[wt_label]
    disp = delta.copy()
    if strip_prefix:
        disp.index = [i.replace(strip_prefix, "") for i in disp.index]
    data = disp.to_numpy(float)
    vmax = np.nanmax(np.abs(data)) if np.isfinite(data).any() else 1.0
    norm = colors.TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)

    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(data, aspect="auto", cmap=cmap, norm=norm)
    ax.set_yticks(np.arange(disp.shape[0])); ax.set_yticklabels(disp.index)
    ax.set_xticks(np.arange(disp.shape[1])); ax.set_xticklabels(disp.columns, rotation=0)
    ax.set_xlabel("Compartment"); ax.set_ylabel("Cell type")
    ax.set_title(title or f"Δ proportion ({mut_label} − {wt_label}) (section-equalized)")
    for spine in ['top', 'bottom', 'left', 'right']:
        ax.spines[spine].set_linewidth(2)
    ax.tick_params(axis='both', which='major', width=2.5, length=6)
    ax.tick_params(axis='both', which='minor', width=0.75, length=3)
    cbar = fig.colorbar(im, ax=ax, pad=0.02)
    cbar.ax.tick_params(length=6, width=3)
    cbar.set_label("Δ posterior abundance (mean, section-equalized)")
    plt.tight_layout()
    return delta, fig, ax

# dependent on pieplot proportions
def delta_lollipop_from_cellprop(
    cell_prop: dict,                 # {'wt1': F_drawn_df, 'wt2': F_drawn_df, 'mgb2': F_drawn_df}
    slide_to_geno: dict,             # {"wt1":"WT","wt2":"WT","mgb2":"mutant"}
    *,
    exclude_slides=None,
    wt_label="WT",
    mut_label="mutant",
    prefix_to_remove="c2l_",
    title="Δ proportion (mutant − WT)",
    figsize=(7,6),
    cmap_name="bwr",                 # used only if base_colors not supplied or for fallbacks
    vmin=None, vmax=None,
    base_colors: dict=None,          # NEW: per–cell type palette, e.g. {"Fib2":"lightblue", "detSMC":"lime", ...}
    default_color="gray"             # used if a type not in base_colors and no cmap fallback
):
    """
    Aggregates per-slide compositions (equal weight per slide), computes delta (mutant − WT),
    and plots a lollipop. If `base_colors` is provided, each stick/marker is colored by that palette.
    Otherwise, a centered blue→white→red colormap is used.
    """
    import matplotlib

    # 0) optional exclusions
    if exclude_slides:
        cell_prop = {k: v for k, v in cell_prop.items() if k not in set(exclude_slides)}

    # 1) per-slide composition (mean over spots; rows sum to 1 across selected cols)
    slide_comp = {}
    for slide, F_drawn in cell_prop.items():
        if not isinstance(F_drawn, pd.DataFrame):
            F_drawn = pd.DataFrame(F_drawn)
        slide_comp[slide] = F_drawn.mean(axis=0)
    slide_comp = pd.DataFrame(slide_comp).T  # rows=slides, cols=cell types

    # 2) map slides → genotype (equal weight per slide)
    try:
        slide_comp["genotype"] = [slide_to_geno[s] for s in slide_comp.index]
    except KeyError as e:
        raise KeyError(f"Slide '{e.args[0]}' missing in slide_to_geno mapping.") from e

    present = set(slide_comp["genotype"].unique())
    if not {wt_label, mut_label}.issubset(present):
        raise ValueError(f"Need both genotypes '{wt_label}' and '{mut_label}' present. Found: {present}")

    by_geno_mean = slide_comp.groupby("genotype").mean(numeric_only=True)
    delta = (by_geno_mean.loc[mut_label] - by_geno_mean.loc[wt_label])

    # ---------- plotting prep ----------
    # Clean labels for display
    vals = delta.copy()
    if prefix_to_remove:
        clean_names = [s.replace(prefix_to_remove, "") for s in vals.index]
    else:
        clean_names = list(vals.index)
    vals.index = clean_names
    vals = vals.sort_values()
    y = np.arange(len(vals))

    fig, ax = plt.subplots(figsize=figsize)

    if base_colors:  # --- categorical coloring from palette ---
        # Build color list per row using clean names
        color_list = []
        for ct in vals.index:
            if ct in base_colors:
                color_list.append(base_colors[ct])
            else:
                # try to map back to original key with prefix if user palette uses original names
                orig = (prefix_to_remove + ct) if prefix_to_remove else ct
                color_list.append(base_colors.get(orig, default_color))
        # Draw
        for yi, v, c in zip(y, vals.values, color_list):
            ax.hlines(yi, 0, v, color=c, linewidth=3)
            ax.plot(v, yi, "o", color=c, markersize=6)
        show_colorbar = False

    else:  # --- diverging cmap centered at 0 (fallback) ---
        absmax = float(np.nanmax(np.abs(vals.values))) if len(vals) else 1.0
        cmap = matplotlib.colormaps[cmap_name]
        if vmin is None and vmax is None:
            norm = colors.TwoSlopeNorm(vmin=-absmax, vcenter=0.0, vmax=absmax)
        else:
            norm = colors.TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)
        colors_arr = cmap(norm(vals.values))
        for yi, v, c in zip(y, vals.values, colors_arr):
            ax.hlines(yi, 0, v, color=c, linewidth=2.5)
            ax.plot(v, yi, "o", color=c, markersize=6)
        # colorbar
        sm = cm.ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, pad=0.02)
        cbar.set_label("Δ proportion")
        show_colorbar = True

    # axes cosmetics
    ax.axvline(0, linestyle="--", linewidth=1, color="gray")
    ax.set_yticks(y); ax.set_yticklabels(vals.index)
    ax.set_xlabel(f"Δ proportion ({mut_label} − {wt_label})")
    ax.set_title(title)
    ax.set_ylim(-0.5, len(vals)-0.5)
    ax.grid(False)
    for spine in ['top', 'bottom', 'left', 'right']:
        ax.spines[spine].set_linewidth(2)
    ax.tick_params(axis='both', which='major', width=2.5, length=6)
    ax.tick_params(axis='both', which='minor', width=0.75, length=3)
    plt.tight_layout()

    delta.name = f"delta_{mut_label}_minus_{wt_label}"
    return delta.sort_values(), fig, ax

def plot_delta_bars(delta, base_colors=None, prefix_to_remove="c2l_", figsize=(6,5),
                    title="Δ proportion (mutant − WT)", sort_by_abs=True):
    vals = delta.copy()
    if prefix_to_remove:
        vals.index = [s.replace(prefix_to_remove, "") for s in vals.index]

    if sort_by_abs:
        vals = vals.reindex(vals.abs().sort_values(ascending=True).index)

    colors = None
    if base_colors:
        colors = [base_colors.get(k, base_colors.get(prefix_to_remove + k, "gray"))
                  for k in vals.index]

    fig, ax = plt.subplots(figsize=figsize)
    ax.barh(vals.index, vals.values, color=colors or "gray")
    ax.axvline(0, color="black", lw=1)
    ax.set_xlabel("Δ proportion (mutant − WT)")
    ax.set_title(title)
    plt.tight_layout()
    return fig, ax

def plot_compartment_delta_heatmap(
    cell_prop: dict,            # {'wt1': F_drawn_df, 'wt2': F_drawn_df, 'mgb2': F_drawn_df}
    meta_data: pd.DataFrame,    # columns: sample_id, barcode_norm, region
    slide_to_geno: dict,        # {'wt1':'WT','wt2':'WT','mgb2':'mutant'}
    *,
    sample_id_col="sample_id",
    barcode_col="barcode_norm",
    region_col="region",
    compartments_order=None,    # e.g. ["Inner","Middle","Outer","Umblical"]
    wt_label="WT",
    mut_label="mutant",
    prefix_to_remove="c2l_",    # strip for display-only
    figsize=(8,6)
):
    """
    Builds a matrix of Δ proportion (mutant − WT) per compartment and plots a heatmap.
    Assumes each F_drawn in cell_prop is (spots x selected_celltypes) with rows summing to 1 across selected cols.
    Returns (delta_df, fig, ax) where delta_df rows=cell types, cols=compartments.
    """

    # --- helper: find best separator that matches F_drawn indices to meta composite keys
    def _best_sep(F_idx: pd.Index, meta_sub: pd.DataFrame):
        # try common seps
        seps = ("_","-","","|",":")
        fset = set(F_idx.astype(str))
        best, hits = None, -1
        for s in seps:
            comp = (meta_sub[sample_id_col].astype(str) + s + meta_sub[barcode_col].astype(str)).values
            h = sum((c in fset) for c in comp)
            if h > hits:
                best, hits = s, h
        if hits == 0:
            raise ValueError("Could not match any meta rows to spot indices; check sample_id/barcode_norm format.")
        return best

    # --- 1) per-slide, per-compartment compositions (mean of spot fractions within each region)
    per_slide_region = {}
    for slide, F in cell_prop.items():
        F = F.copy() if isinstance(F, pd.DataFrame) else pd.DataFrame(F)
        # subset meta_data to this slide's sample_id (the sample_id string contained in obs_names)
        # robust approach: meta_data has explicit sample_id; use that to subset
        # guess this slide's sample_id: any row of F index should start with sample_id + sep
        # just filter meta_data by rows whose sample_id starts with slide (common: 'wt1-...')
        msub = meta_data[meta_data[sample_id_col].astype(str).str.startswith(slide)].copy()
        if msub.empty:
            # fallback: exact match (if slide already equals sample_id)
            msub = meta_data[meta_data[sample_id_col].astype(str) == slide].copy()
        if msub.empty:
            raise ValueError(f"No meta_data rows found for slide '{slide}' in column '{sample_id_col}'.")

        sep = _best_sep(F.index, msub)
        msub["__key"] = msub[sample_id_col].astype(str) + sep + msub[barcode_col].astype(str)
        lut = msub.set_index("__key")[region_col].astype(str)

        # map each spot to a region
        r = pd.Series(F.index.astype(str)).map(lut).values
        mask = pd.notnull(r)
        if not mask.any():
            raise ValueError(f"Could not map regions for slide '{slide}'. Check composite key format.")
        F_use = F.loc[mask].copy()
        regions = pd.Index(r[mask], name="region")

        reg_comp = (
            F_use.assign(region=regions.values)
                .groupby("region")
                .mean(numeric_only=True)
        )
        per_slide_region[slide] = reg_comp

    # --- 2) aggregate per genotype to get means, then Δ = mutant − WT per compartment
    # stack into long table: slide, region, celltype, value
    long_rows = []
    for slide, df in per_slide_region.items():
        df2 = df.copy()
        df2["slide"] = slide
        df2["genotype"] = slide_to_geno[slide]
        long_rows.append(df2.reset_index())
    long_df = pd.concat(long_rows, ignore_index=True)

    cell_cols = [c for c in long_df.columns if c not in {"region","slide","genotype"}]
    regions = sorted(long_df["region"].unique()) if compartments_order is None else compartments_order

    # Build delta matrix
    delta_parts = []
    for comp in regions:
        sub = long_df[long_df["region"] == comp]
        if sub.empty:
            continue
        mean_by_geno = sub.groupby("genotype")[cell_cols].mean(numeric_only=True)
        if (wt_label in mean_by_geno.index) and (mut_label in mean_by_geno.index):
            delta_col = (mean_by_geno.loc[mut_label] - mean_by_geno.loc[wt_label])
            delta_col.name = comp
            delta_parts.append(delta_col)
        else:
            # one genotype missing -> fill NaNs for that column
            delta_parts.append(pd.Series(index=cell_cols, dtype=float, name=comp))

    delta_df = pd.concat(delta_parts, axis=1)
    # pretty row labels (display only)
    if prefix_to_remove:
        delta_df.index = [i.replace(prefix_to_remove, "") for i in delta_df.index]

    # --- 3) plot heatmap with matplotlib (center at 0) ---
    fig, ax = plt.subplots(figsize=figsize)
    data = delta_df.to_numpy(dtype=float)
    # symmetric color around 0 for fair visual comparison
    vmax = np.nanmax(np.abs(data)) if np.isfinite(data).any() else 1.0
    norm = colors.TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    im = ax.imshow(data, aspect="auto", cmap="bwr", norm=norm)

    # ticks & labels
    ax.set_yticks(np.arange(delta_df.shape[0]))
    ax.set_yticklabels(delta_df.index)
    ax.set_xticks(np.arange(delta_df.shape[1]))
    ax.set_xticklabels(delta_df.columns, rotation=0)
    ax.set_xlabel("Compartment")
    ax.set_title("Δ proportion (mutant − WT) per compartment")

    # colorbar
    cbar = fig.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label("Δ proportion")

    plt.tight_layout()
    return delta_df, fig, ax
