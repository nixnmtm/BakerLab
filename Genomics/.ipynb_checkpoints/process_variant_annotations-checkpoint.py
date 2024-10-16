import pandas as pd
import seaborn as sns
import os.path as path
import matplotlib.pyplot as plt
import subprocess
import logging
import numpy as np


def process_cmd(cmd):
    output = subprocess.check_output(cmd.split(), 
                                     stderr=subprocess.DEVNULL, 
                                     universal_newlines=True) # universal_newlines converts the object to str
    return output.strip().split("\n")

def split_GT(df, sample_cols):
    for col in sample_cols:
        a1 = []
        a2 = []
        for i in df[col].values:
            a1.append(i[0])
            a2.append(i[2])
        if len(df[col].values) == len(a1):
            df[f"{col}_A1"] = a1
        else:
            print("Length mismatch")
        if len(df[col].values) == len(a2):    
            df[f"{col}_A2"] = a2
        else:
            print("Length mismatch")
    return df

def postprocess_myanno(file_full_path, sample_id=None):
    """
    file_full_path: (str) absolute path of the annovar multianno.txt file
    sample_id: input the sample id if single sample vcf is used/if multisample vcf
               use sample_id=None
    
    Returns:
    Preprocessed with samples names added and Otherinfo cols fixed
    """
    
    myanno = pd.read_csv(file_full_path, delimiter="\t", index_col="Chr", low_memory=False)
    myanno[['gnomad40_exome_AF']] = myanno[['gnomad40_exome_AF']].apply(pd.to_numeric, errors='coerce')
    myanno.rename(columns={"Otherinfo1": "AFvcf1",
                        "Otherinfo2": "QUALvcf1",
                        "Otherinfo4": "Chromosome",
                       "Otherinfo5": "Startvcf1",
                       "Otherinfo6": "IDvcf1",
                       "Otherinfo7": "REFvcf1",
                       "Otherinfo8": "ATLvcf1",
                       "Otherinfo9": "QUAL",
                       "Otherinfo10": "FILTER",
                        "Otherinfo11": "INFO",
                        "Otherinfo12": "FORMAT"}, errors="raise", inplace=True)
    if sample_id == None: 
        sampl_cols = myanno.loc[:, 'Otherinfo13':].columns
        # get sample names
        try:
            vcfpath = path.splitext(file_full_path)[0]+'.vcf'
            cmd = f"/Users/nxr042/research/setups/bcftools1.19/bcftools query -l {vcfpath}"
            sample_names = process_cmd(cmd=cmd)
        except subprocess.CalledProcessError:
            print("Check the VCF path, please include *multianno.vcf in the same path as *multianno.txt")
        if len(sample_names) > 1:
            if len(sampl_cols) == len(sample_names):
                for i in range(len(sample_names)):
                    myanno.columns = myanno.columns.str.replace(sampl_cols[i], sample_names[i])
            df = split_GT(myanno, sample_cols=sample_names)
        else:
            raise Exception("If single sample VCF is used, please input the sample_id")
    else:
        df = myanno.rename(columns={"Otherinfo13": sample_id}, errors="raise")
        sample_names = df.loc[:, sample_id:].columns.values
    return df, sample_names
    

def get_columns(annovar_dbs_list, all_annovar_cols):
    """
    get the database names and find columns related to it in annovar annotation
    
    annovar_dbs_list : list of databases you are interested in ex: functional score databases
    all_annovar_cols : all columns in a annovar annotation file (annovar.hg38_multianno.txt)
    
    return columns that are corresponding to the database list given
    """
    base_cols = ["Start", "End", "Ref", "Alt", "Gene.refGene", "ExonicFunc.refGene"]
    dbs = [i.lower() for i in annovar_dbs_list]
    _cols = []
    for col in all_annovar_cols:  
        check_ = [True if db in col.lower() else False for db in dbs ]
        if any(check_):
            _cols.append(col)
    return base_cols+_cols

def filter_variants(df, af_cut="0.001", all_type_variants=False):
    """
    df: dataframe containing all the variants
    af_cut: allele frequency cut
    
    """
    from natsort import index_natsorted

    if all_type_variants:
        variants = df[(df["gnomad40_exome_AF"] < float(af_cut))]
    else:
        non_splice_vaiants = df[(df["gnomad40_exome_AF"] < float(af_cut)) & 
            ((df["ExonicFunc.refGene"] != "synonymous SNV") & 
             (df["ExonicFunc.refGene"] != "unknown") &
             (df["ExonicFunc.refGene"] != "."))]
        splice_variants = df[(df["gnomad40_exome_AF"] < float(af_cut)) &
            ((df["Func.refGene"] == "splicing") | 
            (df["Func.refGene"] == "exonic;splicing") | 
            (df["Func.refGene"] == "ncRNA_exonic;splicing") | 
            (df["Func.refGene"] == "ncRNA_splicing"))]
    
        variants = pd.concat([non_splice_vaiants, splice_variants]).reset_index()
        variants.sort_values(by="Chr", key=lambda x: np.argsort(index_natsorted(variants["Chr"])))
        variants.set_index("Chr", inplace=True)
    #causal_variants = causal_variants.sort_index(ascending=True, kind='mergesort')
    return variants
    
def check_smc_genes(df):
    import numpy as np
    gs = df["Gene.refGene"].drop_duplicates().tolist()
    smc_genes = np.intersect1d(gs, smc_gene_list)
    if len(smc_genes) > 0 :
        print(True)
        return smc_genes
    else:
        return False  
    
def get_atleast2deleterious(df):
    """
    CHECK THE TRUTH AND RETURN A DF WITH NEW COLUMN
    
    """
    # replace periods with 0, to set logical operators for the comparison
    df = df.replace('.',0)
    df['CADD_phred'] = df['CADD_phred'].astype(float)
    
    mask = (((df["SIFT_pred"] == "D") & (df["CADD_phred"] >= 20)) | 
            ((df["SIFT_pred"] == "D") & (df["Polyphen2_HDIV_pred"] == "D")) |
            ((df["SIFT_pred"] == "D") & (df["Polyphen2_HDIV_pred"] == "P")) |
            ((df["SIFT_pred"] == "D") & (df["MutationTaster_pred"] == "D")) |
            ((df["SIFT_pred"] == "D") & (df["MutationTaster_pred"] == "A")) |
            ((df["CADD_phred"] >= 20) & (df["Polyphen2_HDIV_pred"] == "D")) |
            ((df["CADD_phred"] >= 20) & (df["Polyphen2_HDIV_pred"] == "P")) | 
            ((df["CADD_phred"] >= 20) & (df["MutationTaster_pred"] == "D")) | 
            ((df["Polyphen2_HDIV_pred"] == "D") & (df["MutationTaster_pred"] == "D")) |
            ((df["Polyphen2_HDIV_pred"] == "P") & (df["MutationTaster_pred"] == "D")) |
            ((df["Polyphen2_HDIV_pred"] == "D") & (df["MutationTaster_pred"] == "A")) |
            ((df["Polyphen2_HDIV_pred"] == "P") & (df["MutationTaster_pred"] == "A")))
    df["atleast_2_deleterious"] = mask
    cols = df.columns
    _index = df.columns.get_loc("SIFT_pred") 
    df.insert(_index, 'atleast_2_deleterious', df.pop('atleast_2_deleterious'))
    return df

#import mygene
def get_geneSym(gene_name_or_id, organism="Homo_sapiens"):
    """
    Parsing id and name from the tsv file obtained from the result of below commands
    
    python /Users/nxr042/research/exome/query_ensemble.py --organism homo_sapiens
    python /Users/nxr042/research/exome/query_ensemble.py --organism mus_musculus
    
    gene_name_or_id: given Ensembl_gene_ID returns the HGNC_symbol and vice versa 
    organism: "homo_sapiens or mus_musculus"
    """
    organism = organism.lower()
    data = pd.read_csv(f"/Users/nxr042/research/exome/{organism}_gene_annotation.tsv", sep="\t")
    ens2sym = dict(zip(data.Ensembl_gene_ID, data.HGNC_symbol))
    sym2ens = dict(zip(data.HGNC_symbol, data.Ensembl_gene_ID))
    if gene_name_or_id.startswith('ENSG00') or gene_name_or_id.startswith('ENSMUSG00'):
        try:
            return ens2sym[gene_name_or_id]
        except:
            #print(f"Ensemble id '{gene_name_or_id}' not found")
            return None
    else:
        try:
            return sym2ens[gene_name_or_id]
        except:
            return None
            #print(f"Gene symbol '{gene_name_or_id}' not found")

def fetch_human_protein_atlas(gene_name):
    import requests as req
    import json

    # retreive information from ensemble

    ens_id = get_geneSym(gene_name, organism="Homo_sapiens")
    #print(f"From atlas fn: {ens_id}")
    url = f'https://www.proteinatlas.org/{ens_id}.json'
    resp = req.get(url)
    if resp.ok:
        gene_details = json.loads(resp.text)
        return gene_details

def fetch_mygene(gene_name):
    import mygene
    # retrive info from mygene
    ens_id = get_geneSym(gene_name, organism="Homo_sapiens")
    #print(f"From mygene fn: {ens_id}")
    mg = mygene.MyGeneInfo()
    g = mg.getgene(ens_id)
    return g

import functools

def log_function(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as exc:
            print(args, kwargs, repr(exc))

    return wrapper

@log_function
def fetch_gene_details(gene_name):
    """
    Retrieve gene details from Human Protein Atlas and MyGene
    gene_name: Ensembl ID or gene name
    """
    d = dict()
    d["Gene.refGene"] = gene_name
    d["Center"] = "NCH"
    
    # Fetch protein atlas details
    protein_atlas = fetch_human_protein_atlas(gene_name)
    if protein_atlas:
        for key in ["RNA tissue cell type enrichment", "Single cell expression cluster", "Tissue expression cluster"]:
            if key in protein_atlas:
                if isinstance(protein_atlas[key], list) and len(protein_atlas[key]) > 1:
                    d[key] = ";".join(protein_atlas[key])
                else:
                    d[key] = protein_atlas[key]
            else:
                d[key] = None
    
    # Fetch gene-specific details from MyGene
    mygene = fetch_mygene(gene_name)
    if mygene:
        d["Summary"] = mygene.get("summary", None)
    
    return d


def populate_gene_details(sample_id, absolute_annovar_path, af_cut=0.001, basic_cols=True, all_variants=True):
    """
    **Populates gene deatils per variant**
    Runs the postprocess in annovar annotated results, filters variants for given allele frequency
    and fetch gene details from human protein atlas and my gene to add along with the annovar results.
    
    Parameters:
    :param: sample_id: patient sample id.
    :param: absolute_annovar_path: path of the annotated file
    :param: af_cut: allele frequency cutoff (default 0.1%)
    :return: a dataframe with basic columns used for interpretation and analysis
    
    # https://superfastpython.com/threadpoolexecutor-concurrent-list-comprehension/
    # Based on the above linked used concurrent 
    
    # In the case of PBS, we have to rely on PolyPhen2 HDIV predcition than PolyPhen2 HVAR --> http://genetics.bwh.harvard.edu/pph2/dokuwiki/overview
    
    """
    
    from concurrent.futures import ThreadPoolExecutor, as_completed, ProcessPoolExecutor
    import time
    
    # fetch annovar annotated results and postprocess it
    vc, sample_cols = postprocess_myanno(absolute_annovar_path, sample_id=sample_id)
    
    # causal variants below aaf cutoff
    variants = filter_variants(vc, str(af_cut), all_type_variants=all_variants)
    
    if basic_cols:
        basic_cols = ["Start", "End","Ref", "ATLvcf1", "Gene.refGene", "ExonicFunc.refGene", "AAChange.refGene", "avsnp150", "gnomad40_exome_AF", 
              "Interpro_domain", "SIFT_pred","Polyphen2_HVAR_pred","Polyphen2_HDIV_pred", "ClinPred_pred", "CADD_phred", "MutationTaster_pred",
              'CLNSIG', 'phastCons100way_vertebrate', 'phastCons30way_mammalian', "FORMAT", sample_id]
        variants = variants[basic_cols].reset_index()
    else:
        variants = variants.reset_index()
    
    # get gene details and merge with causal variants dataframe
    genelist = np.unique(variants["Gene.refGene"].values)
    with ThreadPoolExecutor(4) as tpe:
        # execute tasks in parallel and gather the results
        results_list = [result for result in tpe.map(fetch_gene_details, genelist)]       
    gene_details = pd.DataFrame(results_list)
    
    # check for missing rows
    try:
        mask = gene_details.isnull().any(axis=1)
        nanrows = gene_details[mask].shape[0]
        not_nanrows = gene_details[~mask].shape[0]
        if len(results_list) == (not_nanrows + nanrows):
            print(f"All records processed for {sample_id}")
    except:
        print(f"{len(results_list)-(not_nanrows + nanrows)} records missing")
    else:
        # Merge mygene and proteinatlas data with postprocessed annovar results
        merged = pd.merge(variants,gene_details, on="Gene.refGene", how="outer")
        merged.insert(0, "sample_id", sample_id)
       
        # insert a column checking for atleast 2 deleterious predictions among 
        # SIFT_pred, CADD_phred, Polyphen2_HDIV_pred, MutationTaster_pred, ClinPred_pred
        merged = get_atleast2deleterious(merged)
        
        outpath = path.dirname(absolute_annovar_path)
        if all_variants and base_cols:
            merged.to_csv(path.join(outpath, f"af_{af_cut}_all_variants_imp_annotations.csv"), index=False)
        elif all_variants:
            merged.to_csv(path.join(outpath, f"af_{af_cut}_all_variants_all_annotations.csv"), index=False)
        elif basic_cols and not all_variants:
            merged.to_csv(path.join(outpath, f"af_{af_cut}_exonic_variants_imp_annotations.csv"), index=False)
            
        return merged

def split_info(info):
    """
    Splits a semicolon-separated string into a dictionary.
    
    The input string is expected to contain key-value pairs separated by '='. 
    Each key-value pair is separated by a semicolon. 

    param: info: A semicolon-separated string containing key-value pairs (e.g., "key1=value1;key2=value2").
    return: A dictionary where keys and values are derived from the input string.
    """
    return dict(item.split('=') for item in info.split(';'))

def process_vep_annotations(sample_id, absolute_path=None, write=True):
    """
    Processes VEP (Variant Effect Predictor) annotations for a given sample.

    This function reads a VEP output file, splits the 'Extra' column into multiple 
    columns based on key-value pairs, and returns a new DataFrame containing the 
    original data along with the newly created columns.

    param: sample_id: The ID of the sample being processed.
    param: absolute_path: The absolute file path to the VEP output file (must be a 
                          whitespace-delimited text file).

    return: A DataFrame containing the original VEP data along with additional columns 
            extracted from the 'Extra' column.
    """
    
    vep_data = pd.read_csv(absolute_path, delim_whitespace=True, skiprows=92)
    split_extra = vep_data['Extra'].apply(split_info).apply(pd.Series)
    vep_data = vep_processed_data.drop('Extra', axis=1)
    processed_vep = pd.concat([vep_data, split_extra], axis=1)
    processed_vep["sample_id"] = sample_id
    
    if write:
        outpath = path.dirname(absolute_path)
        filepath = path.join(outpath, f"processed_vep_table.csv")
        processed_vep.to_csv(filepath, index=False)
        print(f"Writing vep processed table to {filepath}")
        
    return processed_vep