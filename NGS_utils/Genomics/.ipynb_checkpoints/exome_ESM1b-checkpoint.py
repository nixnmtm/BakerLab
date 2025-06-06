import re
import pandas as pd
import numpy as np
from zipfile import ZipFile

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import gridspec
import seaborn as sns
import os.path as path
import scipy.stats as stats
from sklearn.metrics import roc_curve
import logging

import requests
import json



"""
# esm-variants for missense
https://huggingface.co/spaces/ntranoslab/esm_variants/blob/main/app.py
# for indels use the tutorial below
https://github.com/ntranoslab/esm-variants/blob/main/esm_variants_utils.ipynb
""" 

# for logging
def setup_logger(name="variant_logger", level="INFO"):
    logger = logging.getLogger(name)
    if not logger.handlers:
        level = getattr(logging, level.upper(), logging.INFO)
        logger.setLevel(level)
        ch = logging.StreamHandler()
        ch.setLevel(level)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        ch.setFormatter(formatter)
        logger.addHandler(ch)
    return logger

logger = setup_logger(level="INFO")

# This function is modified from 
# https://huggingface.co/spaces/ntranoslab/esm_variants/blob/main/app.py
def load_LLR(uniprot_id, llr_file):
    """
    Load Log-Likelihood Ratios (LLRs) for a specific UniProt protein ID.

    Parameters:
    ----------
    uniprot_id : str
        The UniProt ID of the protein (e.g., 'P01116' or 'P01116-2').

    llr_file : str
        The absolute path to the ZIP file containing LLR CSVs.

    Returns:
    -------
    pd.DataFrame
        A DataFrame of shape (20, L), where:
        - Rows correspond to amino acid substitutions in the order:
          ['K', 'R', 'H', 'E', 'D', 'N', 'Q', 'T', 'S', 'C', 
           'G', 'A', 'V', 'L', 'I', 'M', 'P', 'Y', 'F', 'W']
        - Columns are labeled by wild-type amino acid and position (e.g., "G 12").

    Example:
    -------
    >>> load_LLR('P01116', '/path/to/LLRs.zip')
    """
    with ZipFile(llr_file) as myzip:
        name = next((f for f in myzip.namelist() if uniprot_id in f), None)
        if name is None:
            raise ValueError(f"No LLR CSV found for UniProt ID {uniprot_id}")
        with myzip.open(name) as file:
            return pd.read_csv(file, index_col=0)

# This function is taken from 
# https://huggingface.co/spaces/ntranoslab/esm_variants/blob/main/app.py
def meltLLR(LLR,gene_prefix=None,ignore_pos=False):
  vars = LLR.melt(ignore_index=False)
  vars['variant'] = [''.join(i.split(' '))+j for i,j in zip(vars['variable'],vars.index)]
  vars['score'] = vars['value']
  vars = vars.set_index('variant')
  if not ignore_pos:
    vars['pos'] = [int(i[1:-1]) for i in vars.index]
  del vars['variable'],vars['value']
  if gene_prefix is not None:
    vars.index=gene_prefix+'_'+vars.index
  return vars

def fetch_gnomad_variants_by_gene(gene_symbol, dataset="gnomad_r4", reference_genome="GRCh38"):
    """
    Fetches variant data for a gene from gnomAD using GraphQL API.
    """
    query = """
    query ($geneSymbol: String!, $referenceGenome: ReferenceGenomeId!, $dataset: DatasetId!) {
      gene(gene_symbol: $geneSymbol, reference_genome: $referenceGenome) {
        variants(dataset: $dataset) {
          variant_id
          consequence
          hgvsc
          hgvsp
          transcript_consequence {
            transcript_id
            gene_symbol
            consequence_terms
            lof
            lof_filter
            lof_flags
          }
          genome {
            ac
            an
            af
          }
        }
      }
    }
    """

    variables = {
        "geneSymbol": gene_symbol,
        "referenceGenome": reference_genome,
        "dataset": dataset
    }
    
    headers = {"Content-Type": "application/json"}
    response = requests.post(
        "https://gnomad.broadinstitute.org/api",
        headers=headers,
        data=json.dumps({"query": query, "variables": variables}),
    )

    if response.status_code != 200:
        raise Exception("gnomAD query failed: " + response.text)

    variants = response.json()["data"]["gene"]["variants"]
    return pd.DataFrame(variants)
    

def safe_split(variant_id):
    parts = variant_id.split('-')
    if len(parts) == 4:
        return parts
    else:
        return [None]*4
    
def format_gnomad_variants(df):
    # --- Split 'variant_id' into 4 columns ---
    variant_parts = df['variant_id'].apply(safe_split)
    split_df = pd.DataFrame(variant_parts.tolist(), columns=['Chromosome', 'Position', 'Ref', 'Alt'])

    # If lengths don't match, align indexes manually
    split_df.index = df.index  # Ensure correct alignment
    df = pd.concat([df, split_df], axis=1)

    # Convert position to int, skip NaNs
    df['Position'] = pd.to_numeric(df['Position'], errors='coerce')

    # --- Rename & create new columns ---
    df = df.rename(columns={
        'hgvsc': 'Transcript Consequence',
        'hgvsp': 'Protein Consequence',
        'consequence': 'VEP Annotation',
    })

    # Handle genome safely
    df['Allele Frequency'] = df['genome'].apply(
        lambda x: x.get('af') if isinstance(x, dict) else None
    )

    # Add extra gnomAD-style columns
    df['rsIDs'] = None
    df['HGVS Consequence'] = df['Transcript Consequence']

    # Final columns
    columns = [
        'Chromosome', 'Position', 'Ref', 'Alt',
        'rsIDs', 'HGVS Consequence',
        'Protein Consequence', 'Transcript Consequence',
        'VEP Annotation', 'Allele Frequency'
    ]
    
    # Keep only missense variants
    return df[df['VEP Annotation'] == 'missense_variant'][columns]

def compare_gnomad_data(api_df, local_df):

    #print("API df columns:", api_df.columns.tolist())
    #print("Local df columns:", local_df.columns.tolist())

    # Create copies to avoid modifying originals
    api_df = api_df.copy()
    local_df = local_df.copy()

    # Ensure Chromosome is string
    api_df["Chromosome"] = api_df["Chromosome"].astype(str)
    local_df["Chromosome"] = local_df["Chromosome"].astype(str)

    # Normalize Position column
    api_df["Position"] = pd.to_numeric(api_df["Position"], errors='coerce').astype("Int64")
    local_df["Position"] = pd.to_numeric(local_df["Position"], errors='coerce').astype("Int64")

    # Fill missing rsIDs
    api_df["rsIDs"] = api_df["rsIDs"].fillna("").astype(str)
    local_df["rsIDs"] = local_df["rsIDs"].fillna("").astype(str)

    # Merge to identify differences
    merged = pd.merge(
        api_df,
        local_df,
        on=["Chromosome", "Position", "Transcript Consequence", "Protein Consequence"],
        how="outer",
        indicator=True
    )

    diffs = merged[merged["_merge"] != "both"]

    print(f"API gnomAD variants (missense): {len(api_df)}")
    print(f"Local gnomAD variants (missense): {len(local_df)}")
    print(f"Variants unmatched between them: {len(diffs)}")
    return merged

def enrich_api_with_local(api_df, local_df):

    """
    There is some variants not fetched in the local downloaded copy;
    But the Local donwnloaded copy has more data ; meaning it has af and rsIDs for lot of variants that is not fetch in API
    So decided to combine and enrich it
    
    """
    
    # Ensure consistent types
    api_df = api_df.copy()
    local_df = local_df.copy()

    api_df["Chromosome"] = api_df["Chromosome"].astype(str)
    local_df["Chromosome"] = local_df["Chromosome"].astype(str)

    api_df["Position"] = pd.to_numeric(api_df["Position"], errors='coerce').astype("Int64")
    local_df["Position"] = pd.to_numeric(local_df["Position"], errors='coerce').astype("Int64")

    # Merge to bring over rsIDs and Allele Frequency from local
    enriched = pd.merge(
        api_df,
        local_df[["Chromosome", "Position", "Transcript Consequence", "Protein Consequence", "rsIDs", "Allele Frequency"]],
        on=["Chromosome", "Position", "Transcript Consequence", "Protein Consequence"],
        how="left",
        suffixes=("", "_local")
    )

    # Fill missing rsIDs and Allele Frequency from local
    enriched["rsIDs"] = enriched["rsIDs"].combine_first(enriched["rsIDs_local"])
    enriched["Allele Frequency"] = enriched["Allele Frequency"].combine_first(enriched["Allele Frequency_local"])

    # Drop the temp columns
    enriched.drop(columns=["rsIDs_local", "Allele Frequency_local"], inplace=True)

    return enriched

# Function to convert p.Cys24Met format to C24M format
def convert_hgvs_to_variant(hgvs_consequence):
    # Check if the input is a valid string
    if isinstance(hgvs_consequence, str):
        # Regular expression to extract original amino acid, position, and new amino acid
        match = re.match(r'p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})', hgvs_consequence)
        if match:
            # Map of three-letter to one-letter amino acid codes
            aa_dict = {
                'Ala': 'A', 'Cys': 'C', 'Asp': 'D', 'Glu': 'E', 'Phe': 'F',
                'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Lys': 'K', 'Leu': 'L',
                'Met': 'M', 'Asn': 'N', 'Pro': 'P', 'Gln': 'Q', 'Arg': 'R',
                'Ser': 'S', 'Thr': 'T', 'Val': 'V', 'Trp': 'W', 'Tyr': 'Y'
            }
            
            original_aa = match.group(1)
            position = match.group(2)
            new_aa = match.group(3)

            # Convert to the one-letter format (e.g., p.Cys24Met -> C24M)
            variant_format = f"{aa_dict[original_aa]}{position}{aa_dict[new_aa]}"
            return variant_format
        else:
            return None  # Return None if the format doesn't match the expected pattern
    else:
        return None  # Return None if the input is not a valid string

# Function to convert an HGVS column to variant format
def convert_hgvs_column(df, column_name):
    df['variant'] = df[column_name].apply(convert_hgvs_to_variant)
    return df 

def calculate_damaging_fraction(df, score_column, threshold=-5):
    df = df.copy()
    df['is_damaging'] = df[score_column] < threshold
    df['damaging_fraction'] = df[score_column].rank(method='min', ascending=True) / len(df)
    return df

def determine_cutoff_from_af(df, score_col="score", af_col="Allele Frequency", af_threshold=0.01):
    """
    Determine an ESM1b score cutoff that best separates rare from common variants.
    
    Parameters:
    - df: DataFrame containing score and allele frequency
    - score_col: column name of ESM1b scores
    - af_col: column name of Allele Frequency
    - af_threshold: frequency threshold for "rare" (default: 0.01)
    
    Returns:
    - optimal_score_cutoff: score that gives best separation of rare vs common
    """
    df = df[[score_col, af_col]].dropna()
    
    # Define rare = 1 (positive class), common = 0
    y_true = (df[af_col] < af_threshold).astype(int)
    scores = df[score_col]

    # ROC curve: fpr, tpr, thresholds
    fpr, tpr, thresholds = roc_curve(y_true, -scores)  # negate because lower score = more deleterious

    # Find optimal threshold (e.g., Youden's J)
    j_scores = tpr - fpr
    best_idx = np.argmax(j_scores)
    optimal_score_cutoff = thresholds[best_idx]
    
    return -optimal_score_cutoff  # flip back since we negated above

def plot_esm1bVSaf_data(df, protein_name="Random Protein", af_cut=0.001, 
                        annotate_top=5, annotate_variants=None, savepath=None, bins=50):
    """
    Plot ESM1b score vs Allele Frequency scatter plot with histograms.

    Parameters:
    - df: DataFrame containing 'Allele Frequency', 'score', 'damaging_fraction', and 'variant'
    - protein_name: Title and filename prefix
    - af_cut: Frequency threshold (currently unused here)
    - annotate_top: Number of most damaging variants to annotate
    - annotate_variants: List of specific variants to annotate
    - savepath: Folder path to save the figure
    - bins: Number of bins for histograms
    """

    cmap = plt.get_cmap('viridis')

    allele_freq = df["Allele Frequency"].values
    esm1b_scores = df["score"].values
    damaging_scores_fraction = df["damaging_fraction"].values
    variants = df["variant"].values

    fig = plt.figure(figsize=(12, 6))
    gs = gridspec.GridSpec(4, 4, width_ratios=[0.3, 1, 0.05, 0.05], height_ratios=[0.3, 1, 0.05, 0.05])

    # Scatter plot
    ax_main = fig.add_subplot(gs[1:3, 1:3])
    sc = ax_main.scatter(allele_freq, esm1b_scores, c=damaging_scores_fraction, cmap=cmap, norm=Normalize(vmin=0, vmax=1))
    ax_main.set_xscale('log')
    ax_main.set_xlabel('Allele Frequency', fontsize=12)
    ax_main.set_ylabel('ESM1b Scores', fontsize=12)

    # Colorbar
    cbar_ax = fig.add_subplot(gs[1:3, 3])
    cbar = plt.colorbar(sc, cax=cbar_ax)
    cbar.set_label('Fraction of mutations with more damaging scores')

    # Add cutoff line
    cutoff = determine_cutoff_from_af(df[['Allele Frequency', 'score']])
    ax_main.axhline(cutoff, color="red", linestyle="--", label=f"Cutoff = {cutoff:.2f}")

    # Top histogram (Allele Frequency)
    ax_histx = fig.add_subplot(gs[0, 1:3], sharex=ax_main)
    allele_freq_clean = allele_freq[(~np.isnan(allele_freq)) & (allele_freq > 0)]
    if allele_freq_clean.size > 0:
        hist_bins = np.logspace(np.log10(allele_freq_clean.min()), np.log10(allele_freq_clean.max()), bins)
        ax_histx.hist(allele_freq_clean, bins=hist_bins, color='gray', density=True)
        ax_histx.set_xscale('log')
        ax_histx.get_xaxis().set_visible(False)
        ax_histx.set_yticks([])
        ax_histx.set_ylabel('Counts')
        ax_histx.set_title(f'Protein Data: {protein_name}')
    else:
        ax_histx.set_visible(False)

    # Side histogram (ESM1b scores)
    ax_histy = fig.add_subplot(gs[1:3, 0], sharey=ax_main)
    ax_histy.hist(esm1b_scores, bins=bins, orientation='horizontal', color="blue", density=True)
    ax_histy.invert_xaxis()
    ax_histy.set_xlabel('Counts')
    ax_histy.set_ylabel('ESM1b Scores')
    ax_histy.axhline(cutoff, color="red", linestyle="--", label=f"Cutoff = {cutoff:.2f}")



    # Annotate top damaging variants
    bbox = dict(boxstyle='round', facecolor='yellow', edgecolor='black', alpha=0.5)
    top_esm_indices = np.argsort(esm1b_scores)[:annotate_top]
    for idx in top_esm_indices:
        ax_main.annotate(variants[idx], (allele_freq[idx], esm1b_scores[idx]), 
                         textcoords="offset points", xytext=(0, 5), ha='center', fontsize=8, color='red')

    # Annotate specific variants
    if annotate_variants is not None:
        for idx, variant in enumerate(variants):
            if variant in annotate_variants:
                ax_main.annotate(variant, (allele_freq[idx], esm1b_scores[idx]), 
                                 textcoords="offset points", xytext=(0, 5), ha='center',
                                 fontsize=8, color='green', bbox=bbox)

    # Save and show
    if savepath:
        plt.savefig(path.join(savepath, f"{protein_name}_esm1b_plot.png"), dpi=300, bbox_inches='tight', pad_inches=0.2)
    plt.show()
