#!/usr/bin/env python3

# Train cell2location Visium mapping based on NBRegression posteriors obatined from refrence signatures of single cell

import os, argparse, warnings
import numpy as np
import pandas as pd
import anndata as ad
import cell2location
from cell2location.models import Cell2location
import scanpy as sc

def main():
    ap = argparse.ArgumentParser(description="cell2location: load/build model, train, export, save")
    ap.add_argument("--adata-vis", required=True, help="Path to Visium .h5ad")
    ap.add_argument("--inf-aver", required=True, help="Path to signature matrix (genes x celltypes), CSV/TSV; index=SYMBOL")
    ap.add_argument("--run-name", default="cell2loc_run", help="Output/run directory")
    ap.add_argument("--batch-key", default="sample", help="obs column (e.g., condition/sample)")
    ap.add_argument("--feature-col", default="SYMBOL", help="adata_vis.var column with gene symbols")
    ap.add_argument("--n-cells-per-location", type=int, default=10, help="Only used when building model")
    ap.add_argument("--detection-alpha", type=float, default=20.0, help="Only used when building model")
    ap.add_argument("--max-epochs", type=int, default=30000)
    ap.add_argument("--batch-size", default="none", help="'none' for full batch or an int (e.g., 4096)")
    ap.add_argument("--train-size", type=float, default=1.0)
    ap.add_argument("--accelerator", type=str, default="gpu", help="accelerator types (“cpu”, “gpu”, “tpu”, “ipu”, “hpu”, “mps, “auto”)")
    ap.add_argument("--posterior-num-samples", type=int, default=1000)
    ap.add_argument("--posterior-batch-size", default="auto", help="'auto' -> n_obs, or an int")
    args = ap.parse_args()

    run = args.run_name
    os.makedirs(run, exist_ok=True)

    print("Loading Visium:", args.adata_vis)
    adata_vis = ad.read_h5ad(args.adata_vis)
    
    print("Loading posterior inference:", args.inf_aver)
    inf_aver = pd.read_csv(args.inf_aver, index_col=0)

    # Required by scvi-based models in *every* session
    Cell2location.setup_anndata(adata=adata_vis, batch_key=args.batch_key)

    mod = cell2location.models.Cell2location(adata_vis, 
    cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=args.n_cells_per_location,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=args.detection_alpha)


    # Train
    batch_size = None if str(args.batch_size).lower() == "none" else int(args.batch_size)
    print("Training...")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mod.train(
            max_epochs=args.max_epochs,
            batch_size=batch_size,
            train_size=args.train_size,
            accelerator=args.accelerator,
        )

    model_trained_dir = os.path.join(run)
    mod.save(model_trained_dir, overwrite=True)
    print("Saved trained model to:", model_trained_dir)
    
    # If needed to load the model 
    #adata_file = f"{run}/adata_vis_intersected.h5ad"
    #adata_vis = sc.read_h5ad(adata_file)
    #mod = cell2location.models.Cell2location.load(f"{run}", adata_vis)    

    # Export posterior into adata_vis
    pb = adata_vis.n_obs if str(args.posterior_batch_size).lower() == "auto" else int(args.posterior_batch_size)
    print("Exporting posterior...")
    adata_vis = mod.export_posterior(
        adata_vis,
        sample_kwargs={"num_samples": args.posterior_num_samples, "batch_size": pb},
    )

    # Save AnnData with results
    adata_file = os.path.join(run, "sp.h5ad")
    adata_vis.write(adata_file)
    print("Saved posterior-added AnnData to:", adata_file)

if __name__ == "__main__":
    main()

