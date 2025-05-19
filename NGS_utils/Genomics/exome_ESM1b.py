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
        rsid
        consequence
        hgvsc
        hgvsp
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

    print("API df columns:", api_df.columns.tolist())
    print("Local df columns:", local_df.columns.tolist())

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

    print(f"API variants (missense): {len(api_df)}")
    print(f"Local variants (missense): {len(local_df)}")
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

def plot_hist_kde_with_thresholds(ax, data, bins=30, lower_pct=30, upper_pct=70,
                                 hist_color='gray', kde_color='blue',
                                 lower_line_color='red', upper_line_color='green'):
    """
    Plot a horizontal histogram with KDE on ax, add horizontal threshold lines at percentiles.

    Parameters:
    - ax: matplotlib Axes object where plot will be drawn
    - data: 1D array-like numeric data (e.g., esm1b_scores)
    - bins: number of histogram bins
    - lower_pct: percentile for lower threshold line (e.g., 30 for 30th percentile)
    - upper_pct: percentile for upper threshold line (e.g., 70 for 70th percentile)
    - hist_color: color of the histogram bars
    - kde_color: color of the KDE line
    - lower_line_color: color of the lower threshold horizontal line
    - upper_line_color: color of the upper threshold horizontal line
    """

    # Calculate threshold values from percentiles
    lower_thresh = np.percentile(data, lower_pct)
    upper_thresh = np.percentile(data, upper_pct)

    # Plot normalized horizontal histogram
    ax.hist(data, bins=bins, orientation='horizontal', color=hist_color, density=True)
    ax.invert_xaxis()

    # Calculate KDE for the data
    kde = stats.gaussian_kde(data)
    score_range = np.linspace(min(data), max(data), 1000)
    kde_values = kde(score_range)

    # Plot KDE line
    ax.plot(kde_values, score_range, color=kde_color, linewidth=1.5)

    # Plot horizontal threshold lines
    ax.axhline(lower_thresh, color=lower_line_color, linestyle='--', linewidth=1.5,
               label=f'{lower_pct}th percentile ({lower_thresh:.2f})')
    ax.axhline(upper_thresh, color=upper_line_color, linestyle='--', linewidth=1.5,
               label=f'{upper_pct}th percentile ({upper_thresh:.2f})')

    # Optionally add legend
    ax.legend(loc='upper right')

    # Set label for the y-axis (score axis)
    ax.set_ylabel('ESM1b Scores')


# Function to plot the data
def plot_esm1bVSaf_data(df, protein_name="Random Protein", annotate_top=5, annotate_variants=None, savepath=None):
    """
    Plot ESM1B score vs Allele Frequency scatter plot with histograms

    :param df: DataFrame with 
    df: the merged df that has gnomad data, esm1b scores and variants     
    
    """

    cmap = plt.get_cmap('viridis')
    
    allele_freq= df["Allele Frequency"].values
    esm1b_scores = df["score"].values
    damaging_scores_fraction = df["damaging_fraction"].values
    variants = df["variant"].values

    # Create a scatter plot with histograms for both axes
    fig = plt.figure(figsize=(12,6))
    gs = gridspec.GridSpec(4, 4, width_ratios=[0.3, 1, 0.05, 0.05], height_ratios=[0.3, 1, 0.05, 0.05])
    grid = plt.GridSpec(4, 4, hspace=0.1, wspace=0.6)
    
    # Main scatter plot
    ax_main = fig.add_subplot(gs[1:3, 1:3])
    sc = ax_main.scatter(allele_freq, esm1b_scores, c=damaging_scores_fraction, cmap=cmap, norm=Normalize(vmin=0, vmax=1))
    ax_main.set_xlabel('Allele Frequency', fontsize=12)
    ax_main.set_ylabel("ESM1b Scores", fontsize=12)
    ax_main.set_xscale('log')
    
    # Color bar to the side
    cbar_ax = fig.add_subplot(gs[1:3, 3])
    cbar = plt.colorbar(sc, cax=cbar_ax)
    cbar.set_label('Fraction of mutations with more damaging scores')
    
    # Histograms on top and side 
    ax_histx = fig.add_subplot(gs[0, 1:3], sharex=ax_main)
    ax_histy = fig.add_subplot(gs[1:3, 0], sharey=ax_main)
    
    # Top histogram for ESM1b scores
    ax_histx.hist(allele_freq, bins=np.logspace(np.log10(min(allele_freq)), np.log10(max(allele_freq)), 30), color='gray')
    ax_histx.set_xscale('log')  # Ensure that the histogram x-axis matches the log scale of the scatter plot
    ax_histx.set_title(f'Protein Data: {protein_name}')
    ax_histx.get_xaxis().set_visible(False)  # Hide x-axis labels to prevent clutter
    ax_histx.set_yticks([])  # Remove y-axis ticks from the top histogram
    ax_histx.set_ylabel('Counts')
    ax_histx.set_title(f'Protein Data: {protein_name}')
    
    # Left histogram for ESM1b scores
    plot_hist_kde_with_thresholds(ax_histy, esm1b_scores, lower_pct=40, upper_pct=80)
    
    # Label the variants that are in the top 10 ESM1b scores and top 5 allele frequencies
    top_10_esm_indices = np.argsort(esm1b_scores)[:annotate_top]  # Top 10 highest ESM1b scores
    top_5_allele_indices = np.argsort(allele_freq)[-annotate_top:]  # Top 5 highest allele frequencies

    bbox = dict(dict(boxstyle='round',facecolor='yellow', edgecolor='black', alpha=0.5))
    
    # Annotate the top variants
    for idx in top_10_esm_indices:
        ax_main.annotate(variants[idx], (allele_freq[idx], esm1b_scores[idx]), 
                         textcoords="offset points", xytext=(0,5), ha='center', fontsize=8, color='red')#, bbox=bbox)
    
    for idx in top_5_allele_indices:
        ax_main.annotate(variants[idx], (allele_freq[idx], esm1b_scores[idx]), 
                         textcoords="offset points", xytext=(0,5), ha='center', fontsize=8, color='blue')#, bbox=bbox)
        
    if annotate_variants is not None:
        for idx, variant in enumerate(variants):
            if variant in annotate_variants:
                ax_main.annotate(variant, (allele_freq[idx], esm1b_scores[idx]), 
                                 textcoords="offset points", xytext=(0, 5), ha='center', fontsize=8, color='green', bbox=bbox)
    #fig.tight_layout()
    plt.savefig(path.join(savepath, f"{protein_name}_esm1b_plot.png"), dpi=300, bbox_inches='tight', pad_inches=0.2)
    plt.show()