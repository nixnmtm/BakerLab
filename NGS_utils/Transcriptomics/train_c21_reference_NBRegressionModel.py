#!/usr/bin/env python3
"""
End-to-end Cell2location reference pipeline (HPC/GPU-ready).

Estimation of reference cell type signatures (NB regression)

Steps:
 1) Load AnnData
 2) Remove existing embeddings (UMAP/PCA), optional cleanup
 3) Filter genes using cell2location.utils.filtering.filter_genes
 4) Create var['SYMBOL']
 5) Map long celltype names to short labels; drop unwanted types (e.g., Allantois)
 6) Basic QC metrics
 7) Setup & train RegressionModel (labels=celltype, batch=day)
 8) Export posterior + signatures; save outputs

Example:
  python c2l_ref_full_hpc.py \
    --input merged_bladder_candidates.h5ad \
    --run-name outputs_c2l_ref_full \
    --labels-key celltype \
    --batch-key day \
    --remove-celltypes Allantois \
    --remove-embeddings True \
    --cell-count-cutoff 5 \
    --cell-percentage-cutoff2 0.03 \
    --nonz-mean-cutoff 1.12
"""

import os
import sys
import argparse
import logging
import gc

import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import issparse, csr_matrix

import scanpy as sc
import torch
import cell2location as c2l
from cell2location.utils.filtering import filter_genes

import matplotlib.pyplot as plt

# --------- CLI ---------
def parse_args():
    p = argparse.ArgumentParser(description="Cell2location reference pipeline on HPC (GPU-first).")
    p.add_argument("--input", required=True, help="Path to merged_bladder_candidates.h5ad")
    p.add_argument("--run-name", default="outputs_c2l_ref_full", help="Output directory")
    p.add_argument("--labels-key", default="celltype", help="obs column for cell type labels")
    p.add_argument("--batch-key", default="day", help="obs column for batch")
    p.add_argument("--layer", default="counts", help="counts layer name (mirrors .X if missing)")
    p.add_argument("--remove-celltypes", nargs="+", default=[],
                   help="Cell types to drop (space-separated list)")
    p.add_argument("--remove-embeddings",action='store_true', help="remove the embeddings in the adata")
    # Gene filtering thresholds (Cell2location style)
    p.add_argument("--cell-count-cutoff", type=int, default=5,
                   help="Absolute min #cells expressing a gene (for rare but strong genes)")
    p.add_argument("--cell-percentage-cutoff2", type=float, default=0.03,
                   help="Keep genes seen in at least this fraction of cells (e.g., 0.03 = 3%)")
    p.add_argument("--nonz-mean-cutoff", type=float, default=1.12,
                   help="Min mean count in non-zero cells (linear scale)")
    # Training
    p.add_argument("--max-epochs", type=int, default=250)
    p.add_argument("--batch-size", type=int, default=1024)
    p.add_argument("--posterior-samples", type=int, default=1000)
    p.add_argument("--posterior-batch-size", type=int, default=2500)
    # Misc
    p.add_argument("--drop-umap-obs", action="store_true", help="Drop UMAP_1/2/3 columns in .obs")
    return p.parse_args()


# --------- Utils ---------
def logmem(msg, log):
    try:
        import psutil
        rss = psutil.Process().memory_info().rss / 1e9
        log.info(f"{msg} | RSS ~ {rss:.2f} GB")
    except Exception:
        log.info(msg)

def reclaim_mem():
    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()


# --------- Main ---------
def main():
    args = parse_args()
    
    obs_col = args.labels_key    

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
    )
    log = logging.getLogger("c2l_ref_full")


    # --- Load ---
    log.info("Reading: %s", args.input)
    adata = sc.read_h5ad(args.input)
    adata.obs_names_make_unique()
    log.info("Loaded shape: %s cells x %s genes", adata.n_obs, adata.n_vars)

    # --- Remove embeddings (we’ll compute later if needed) ---
    # .obsm keys
    if args.remove_embeddings:
        for k in list(adata.obsm.keys()):
            if k in ("X_umap", "X_pca"):
                del adata.obsm[k]
    # .obs UMAP columns if requested
    if args.drop_umap_obs:
        cols_to_remove = [c for c in ["UMAP_1", "UMAP_2", "UMAP_3"] if c in adata.obs.columns]
        if cols_to_remove:
            adata.obs.drop(columns=cols_to_remove, inplace=True)
            log.info("Dropped obs columns: %s", cols_to_remove)

    # --- Ensure counts layer exists (mirror .X if needed) ---
    if args.layer not in adata.layers:
        log.warning("Layer '%s' not found, mirroring .X into it (assumed raw counts).", args.layer)
        adata.layers[args.layer] = adata.X.copy()
    # Store counts as CSR
    if not issparse(adata.layers[args.layer]):
        adata.layers[args.layer] = csr_matrix(adata.layers[args.layer])

    # remove clusters with cells less than 3
    vc = adata.obs[obs_col].value_counts()
    keep = vc[vc > 3].index
    if len(keep) < len(vc):
        print("Dropping labels with <3 cells:", sorted(set(vc.index) - set(keep)))
        adata = adata[adata.obs[obs_col].isin(keep)].copy()
    
    # --- Gene filtering (Cell2location recommended) ---
    log.info("Filtering genes with cell2location.utils.filtering.filter_genes(...)")
    selected = filter_genes(
        adata,
        cell_count_cutoff=args.cell_count_cutoff,
        cell_percentage_cutoff2=args.cell_percentage_cutoff2,
        nonz_mean_cutoff=args.nonz_mean_cutoff
    )
    
    run_dir = args.run_name
    os.makedirs(run_dir, exist_ok=True)
    
    plt.savefig(os.path.join(run_dir,"gene_filtering_diagnostic.png"), dpi=300, bbox_inches="tight")
    plt.close()

    log.info("Selected genes: %d / %d", int(len(selected)), adata.n_vars)
    adata_ref = adata[:, selected].copy()
    adata_ref.obs_names_make_unique()
    del adata
    reclaim_mem()
    logmem("After gene filtering", log)

    if obs_col not in adata_ref.obs.columns:
        raise ValueError(f"Expected obs[{obs_col}] in the input AnnData.")

    # Remove specified cell types (list)
    if args.remove_celltypes is not None:
        to_remove = [c for c in args.remove_celltypes if c in adata_ref.obs["celltype_long"].unique()]
        if to_remove:
            for c in to_remove:
                n_drop = (adata_ref.obs["celltype_long"] == c).sum()
                log.info("Removing '%s' (%d cells).", c, n_drop)
            adata_ref = adata_ref[~adata_ref.obs["celltype_long"].isin(to_remove)].copy()
    
    # Ensure categorical dtypes for model keys
    adata_ref.obs[obs_col] = adata_ref.obs[obs_col].astype("category")
    if args.batch_key not in adata_ref.obs.columns:
        raise ValueError(f"batch_key '{args.batch_key}' not found in obs.")
    adata_ref.obs[args.batch_key] = adata_ref.obs[args.batch_key].astype("category")

    # --- Setup for RegressionModel ---
    c2l.models.RegressionModel.setup_anndata(
        adata=adata_ref, 
        batch_key=args.batch_key,
        labels_key=obs_col
    )

    # --- Train ---
    log.info("Initializing RegressionModel…")
    mod = c2l.models.RegressionModel(adata_ref)
    
    # Hopper-friendly math + precision env override
    torch.backends.cuda.matmul.allow_tf32 = True
    torch.set_float32_matmul_precision("high")

    accelerator = "gpu" if torch.cuda.is_available() else "cpu"
    precision = os.environ.get("SCVI_PRECISION", "bf16-mixed") if accelerator == "gpu" else "32-true"

    print(
    "CUDA avail:", torch.cuda.is_available(),
    "| device:", torch.cuda.get_device_name(0) if torch.cuda.is_available() else "CPU",
    "| precision:", precision)

    mod.train(
    max_epochs=250,
    batch_size=1024,
    accelerator=accelerator,
    precision=precision)

    # --- Export posterior & signatures ---
    log.info("Exporting posterior (samples=%d, batch_size=%d)…",
             args.posterior_samples, args.posterior_batch_size)
    adata_ref = mod.export_posterior(
        adata_ref,
        sample_kwargs={"num_samples": args.posterior_samples, "batch_size": args.posterior_batch_size}
    )

    if "means_per_cluster_mu_fg" not in adata_ref.varm:
        raise KeyError("Expected 'means_per_cluster_mu_fg' in adata.varm after export_posterior.")
    sig = pd.DataFrame(
        adata_ref.varm["means_per_cluster_mu_fg"],
        index=adata_ref.var_names,
        columns=adata_ref.uns["mod"]["factor_names"]
    )

    # --- Save outputs ---

    sig_path = os.path.join(run_dir, "cell2location_reference_signatures.csv")
    ad_path  = os.path.join(run_dir, "sc.h5ad")

    log.info("Saving signatures: %s", sig_path)
    sig.to_csv(sig_path)

    log.info("Saving trained model to: %s", run_dir)
    mod.save(run_dir, overwrite=True)

    log.info("Saving .h5ad with posterior: %s", ad_path)
    adata_ref.write_h5ad(ad_path, compression="gzip")

    log.info("Done. Signatures shape: %s", sig.shape)


if __name__ == "__main__":
    main()

