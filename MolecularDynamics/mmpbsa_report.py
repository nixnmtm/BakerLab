#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from io import StringIO


# =============================================================================
# Plot style (publication, no grids)
# =============================================================================
def style_axes(ax):
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", which="both", labelsize=11)


# =============================================================================
# Parsing gmx_MMPBSA outputs
# =============================================================================
DELTA_SECTION_RE = re.compile(
    r"Delta \(Complex - Receptor - Ligand\):\s*(.*?)\n\s*-{5,}",
    flags=re.S | re.M,
)

DEFAULT_COMPONENT_ORDER = ["ΔVDWAALS", "ΔEEL", "ΔEPB", "ΔENPOLAR", "ΔGGAS", "ΔGSOLV", "ΔTOTAL"]

def parse_delta_from_dat(dat_path: str | Path) -> pd.DataFrame:
    """
    Parse the Delta (Complex - Receptor - Ligand) table from FINAL_RESULTS_MMPBSA_*.dat.

    Returns a DataFrame with:
      component (e.g. ΔVDWAALS), avg, sd, sem
    """
    dat_path = Path(dat_path)
    text = dat_path.read_text(errors="replace")

    # Capture everything after the Delta header up to the next dashed separator block end
    m = re.search(
        r"Delta \(Complex - Receptor - Ligand\):\s*\n"
        r"Energy Component.*?\n-+\n"
        r"(.*?)\n-+\s*\n-+",
        text,
        flags=re.S,
    )
    if not m:
        raise ValueError(f"Could not find Delta table in: {dat_path}")

    block = m.group(1)

    rows = []
    for line in block.splitlines():
        line = line.strip()
        if not line:
            continue

        parts = re.split(r"\s+", line)
        comp = parts[0]

        # Normalize component label a bit
        # (keep unicode Δ if present; if user ever sees 'DELTA' style, you can expand this)
        if comp.upper() == "DELTA" and len(parts) > 1:
            comp = "Δ" + parts[1]

        # Expect: Average SD(Prop.) SD SEM(Prop.) SEM  (but we only need avg, sd, sem)
        nums = []
        for tok in parts[1:]:
            try:
                nums.append(float(tok))
            except ValueError:
                pass

        if len(nums) >= 1:
            avg = nums[0]
            sd = nums[2] if len(nums) >= 3 else np.nan
            sem = nums[4] if len(nums) >= 5 else np.nan
            rows.append({"component": comp, "avg": avg, "sd": sd, "sem": sem})

    df = pd.DataFrame(rows)
    if df.empty:
        raise ValueError(f"Delta table found but no rows parsed from: {dat_path}")
    return df

def load_energy_csv(csv_path: str | Path, *, table: str = "complex") -> pd.DataFrame:
    """
    Load gmx_MMPBSA per-frame CSV-like outputs that contain leading section text.
    It finds the requested table header and parses from there.

    table: "complex" | "receptor" | "ligand" | "delta"
    """
    csv_path = Path(csv_path)
    text = csv_path.read_text(errors="replace").splitlines()

    # Identify which block you want
    # Common patterns in gmx_MMPBSA outputs:
    #   "Complex Energy Terms" / "Receptor Energy Terms" / "Ligand Energy Terms"
    #   then "Frame #,...."
    table_key = {
        "complex": "Complex Energy Terms",
        "receptor": "Receptor Energy Terms",
        "ligand": "Ligand Energy Terms",
        "delta": "Delta Energy Terms",
    }[table.lower()]

    # Find the block start (table title), then find the next "Frame #," line
    start_idx = None
    for i, line in enumerate(text):
        if line.strip() == table_key:
            start_idx = i
            break
    if start_idx is None:
        # fallback: just find the first Frame header
        for i, line in enumerate(text):
            if line.strip().startswith("Frame #,"):
                start_idx = i - 1
                break
        if start_idx is None:
            raise ValueError(f"Could not find '{table_key}' or any 'Frame #,' header in {csv_path}")

    header_idx = None
    for j in range(start_idx, len(text)):
        if text[j].strip().startswith("Frame #,"):
            header_idx = j
            break
    if header_idx is None:
        raise ValueError(f"Found '{table_key}' but no 'Frame #,' header after it in {csv_path}")

    # Collect CSV lines until the next blank line or next major section header
    # We stop when a line does not contain commas but is non-empty (new section),
    # or when we hit a known section label.
    data_lines = [text[header_idx]]
    for k in range(header_idx + 1, len(text)):
        s = text[k].strip()
        if not s:
            break
        if "," not in s:
            break
        data_lines.append(text[k])

    df = pd.read_csv(StringIO("\n".join(data_lines)))

    # Clean column names (optional)
    df.columns = [c.strip().replace(" ", "_") for c in df.columns]

    # Convert numeric columns
    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="ignore")

    return df

def load_decomp_perframe_csv(csv_path: str | Path) -> pd.DataFrame:
    """
    Parse FINAL_DECOMP_MMPBSA_*.csv (per-frame decomposition) which contains
    leading text headers and may include CRLF artifacts (^M).
    Returns a tidy DataFrame with numeric columns.
    """
    csv_path = Path(csv_path)

    lines = csv_path.read_text(errors="replace").splitlines()

    # Find the real CSV header line
    header_idx = None
    for i, line in enumerate(lines):
        s = line.strip().replace("\r", "")
        if s.startswith("Frame #,") and "Residue" in s:
            header_idx = i
            break

    if header_idx is None:
        raise ValueError(f"Could not find 'Frame #,Residue,...' header in {csv_path}")

    # Keep lines that look like CSV rows (contain commas) until a blank or a new section
    data_lines = [lines[header_idx].replace("\r", "")]
    for j in range(header_idx + 1, len(lines)):
        s = lines[j].strip().replace("\r", "")
        if not s:
            break
        if "," not in s:
            break
        data_lines.append(lines[j].replace("\r", ""))

    df = pd.read_csv(StringIO("\n".join(data_lines)))

    # normalize column names
    df.columns = [c.strip().replace(" ", "_") for c in df.columns]

    # numeric conversion
    for c in df.columns:
        if c not in ("Residue", "Residue_id", "ResidueID"):
            df[c] = pd.to_numeric(df[c], errors="ignore")

    return df

def load_decomp_summary_dat(dat_path: str | Path) -> pd.DataFrame:
    """
    Parse FINAL_DECOMP_MMPBSA_*.dat summary table (Avg/SD/SEM per residue).
    Returns DataFrame with columns like Internal_Avg, Internal_SD, Internal_SEM, ... TOTAL_Avg, TOTAL_SD, TOTAL_SEM.
    """
    dat_path = Path(dat_path)
    lines = dat_path.read_text(errors="replace").splitlines()

    # find header line
    header1_idx = None
    for i, line in enumerate(lines):
        s = line.strip()
        if s.startswith("Residue,Internal"):
            header1_idx = i
            break
    if header1_idx is None:
        raise ValueError(f"Could not find decomposition header in {dat_path}")

    header1 = [x.strip() for x in lines[header1_idx].split(",")]
    header2 = [x.strip() for x in lines[header1_idx + 1].split(",")]

    # Build grouped names:
    # header1 has repeated empty strings due to grouping (Internal,,,van der Waals,,, ...)
    # header2 has Avg, SD, SEM repeated in order.
    groups = []
    current = None
    for h in header1:
        if h:
            current = h
        groups.append(current)

    # Compose final column names
    cols = []
    for g, h2 in zip(groups, header2):
        if g is None:
            cols.append(h2 or "Unknown")
        else:
            if g == "Residue":
                cols.append("Residue")
            else:
                # standardize the statistic label
                stat = h2.replace("Std. Dev.", "SD").replace("Std. Err. of Mean", "SEM").replace("Avg.", "Avg")
                stat = stat.replace(" ", "_").replace(".", "")
                cols.append(f"{g}_{stat}" if stat else g)

    # Data lines start after the second header row
    data = []
    for j in range(header1_idx + 2, len(lines)):
        s = lines[j].strip()
        if not s:
            continue
        if s.startswith(("DELTAS", "Total", "Side", "Complex", "Receptor", "Ligand")):
            continue
        if "," not in s:
            continue
        data.append(lines[j])

    df = pd.read_csv(StringIO("\n".join(data)), header=None, names=cols)

    # Coerce numeric columns except residue
    for c in df.columns:
        if c != "Residue":
            df[c] = pd.to_numeric(df[c], errors="coerce")

    return df

def guess_residue_id_column(df: pd.DataFrame) -> str:
    candidates = ["Residue", "residue", "resname", "res", "Res", "RES", "ID", "id"]
    for c in candidates:
        if c in df.columns:
            return c
    for c in df.columns:
        if not pd.api.types.is_numeric_dtype(df[c]):
            return c
    raise ValueError("Could not infer residue identifier column in decomposition CSV.")


def pick_total_column(df: pd.DataFrame) -> str:
    # gmx_MMPBSA versions differ; try likely names
    candidates = ["TOTAL", "Total", "ΔTOTAL", "dTOTAL", "DELTA_TOTAL"]
    for c in candidates:
        if c in df.columns and pd.api.types.is_numeric_dtype(df[c]):
            return c
    # fallback: last numeric col
    num_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
    if not num_cols:
        raise ValueError("No numeric columns found in decomposition CSV.")
    return num_cols[-1]


# =============================================================================
# Data model + discovery
# =============================================================================
@dataclass(frozen=True)
class RepFiles:
    dat: Path
    csv: Path
    decomp_dat: Path
    decomp_csv: Path


def expected_files(rep_dir: Path) -> RepFiles:
    rep = rep_dir.name  # prod_1 etc
    return RepFiles(
        dat=rep_dir / f"FINAL_RESULTS_MMPBSA_{rep}.dat",
        csv=rep_dir / f"FINAL_RESULTS_MMPBSA_{rep}.csv",
        decomp_dat=rep_dir / f"FINAL_DECOMP_MMPBSA_{rep}.dat",
        decomp_csv=rep_dir / f"FINAL_DECOMP_MMPBSA_{rep}.csv",
    )


def discover_isoforms(root: Path, only: Optional[List[str]] = None) -> List[Path]:
    if only:
        return [root / x for x in only if (root / x).is_dir()]
    return sorted([p for p in root.iterdir() if p.is_dir()])


def discover_systems(iso_dir: Path, only: Optional[List[str]] = None) -> List[Path]:
    if only:
        return [iso_dir / s for s in only if (iso_dir / s).is_dir()]
    return sorted([
        p for p in iso_dir.iterdir()
        if p.is_dir() and p.name not in {"gmx_MMPBSA"} and not p.name.startswith("prod_")
    ])

def discover_replicates(system_dir: Path) -> List[Path]:
    """
    Your layout: <isoform>/<system>/gmx_MMPBSA/prod_*/
    """
    mmpbsa_dir = system_dir / "gmx_MMPBSA"
    if not mmpbsa_dir.is_dir():
        return []
    return sorted([p for p in mmpbsa_dir.iterdir() if p.is_dir() and p.name.startswith("prod_")])



# =============================================================================
# Stats across replicates
# =============================================================================
def summarize_components_across_reps(rep_component_tables: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """
    Input: dict rep_name -> df(component, avg)
    Output: df(component, mean, sd, sem) across replicates.
    """
    big = []
    for rep, df in rep_component_tables.items():
        tmp = df[["component", "avg"]].copy()
        tmp["replicate"] = rep
        big.append(tmp)
    big = pd.concat(big, ignore_index=True)

    piv = big.pivot_table(index="component", columns="replicate", values="avg", aggfunc="first")
    out = pd.DataFrame({
        "component": piv.index,
        "mean": piv.mean(axis=1),
        "sd": piv.std(axis=1, ddof=1),
    }).reset_index(drop=True)
    out["sem"] = out["sd"] / np.sqrt(piv.shape[1])
    return out


def paired_delta_delta_components(
    wt_rep_tables: Dict[str, pd.DataFrame],
    mut_rep_tables: Dict[str, pd.DataFrame],
) -> pd.DataFrame:
    """
    Compute paired ΔΔ(component) = MUT - WT for each matching replicate, then mean±SEM across pairs.
    """
    common_reps = sorted(set(wt_rep_tables.keys()) & set(mut_rep_tables.keys()))
    if not common_reps:
        raise ValueError("No matching replicate names between WT and MUT (e.g., prod_1).")

    per_pair = []
    for rep in common_reps:
        wt = wt_rep_tables[rep].set_index("component")["avg"]
        mut = mut_rep_tables[rep].set_index("component")["avg"]
        common_components = wt.index.intersection(mut.index)

        dd = (mut.loc[common_components] - wt.loc[common_components]).rename("dd")
        tmp = dd.reset_index().rename(columns={"index": "component"})
        tmp["replicate"] = rep
        per_pair.append(tmp)

    per_pair = pd.concat(per_pair, ignore_index=True)

    piv = per_pair.pivot_table(index="component", columns="replicate", values="dd", aggfunc="first")
    out = pd.DataFrame({
        "component": piv.index,
        "mean_dd": piv.mean(axis=1),
        "sd_dd": piv.std(axis=1, ddof=1),
    }).reset_index(drop=True)
    out["sem_dd"] = out["sd_dd"] / np.sqrt(piv.shape[1])
    out = out.rename(columns={"component": "component"})
    return out


# =============================================================================
# Plotting
# =============================================================================
def plot_components_bar(df: pd.DataFrame, outpng: Path, title: str, value_col="mean", err_col="sem",
                        order: Optional[List[str]] = None, ylabel="Energy (kcal/mol)"):
    order = order or DEFAULT_COMPONENT_ORDER
    dfi = df.set_index("component")
    keep = [c for c in order if c in dfi.index]

    vals = dfi.loc[keep, value_col].values
    errs = dfi.loc[keep, err_col].values if err_col in dfi.columns else None

    fig, ax = plt.subplots(figsize=(7.6, 4.0))
    x = np.arange(len(keep))
    ax.bar(x, vals, yerr=errs, capsize=4 if errs is not None else 0)
    ax.set_xticks(x)
    ax.set_xticklabels(keep, rotation=45, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    style_axes(ax)
    plt.tight_layout()
    fig.savefig(outpng, dpi=300)
    plt.close(fig)


def plot_total_only(df: pd.DataFrame, outpng: Path, title: str, value_col="mean", err_col="sem"):
    dfi = df.set_index("component")
    if "ΔTOTAL" not in dfi.index:
        raise ValueError("ΔTOTAL not found for plot_total_only.")

    mean = float(dfi.loc["ΔTOTAL", value_col])
    err = float(dfi.loc["ΔTOTAL", err_col]) if err_col in dfi.columns and pd.notna(dfi.loc["ΔTOTAL", err_col]) else None

    fig, ax = plt.subplots(figsize=(4.0, 4.0))
    ax.bar([0], [mean], yerr=[err] if err is not None else None, capsize=4 if err is not None else 0)
    ax.set_xticks([0])
    ax.set_xticklabels(["ΔTOTAL"])
    ax.set_ylabel("Energy (kcal/mol)")
    ax.set_title(title)
    style_axes(ax)
    plt.tight_layout()
    fig.savefig(outpng, dpi=300)
    plt.close(fig)


def plot_timeseries(csv_path: Path, outdir: Path, prefix: str,
                    cols=("TOTAL", "GGAS", "GSOLV", "EEL", "VDWAALS", "EPB", "ENPOLAR")):
    df = load_energy_csv(csv_path)
    for col in cols:
        if col not in df.columns:
            continue
        fig, ax = plt.subplots(figsize=(7.6, 3.2))
        ax.plot(df[col].values, linewidth=1)
        ax.set_xlabel("Frame index")
        ax.set_ylabel(f"{col} (kcal/mol)")
        ax.set_title(f"{prefix}: {col} vs frame")
        style_axes(ax)
        plt.tight_layout()
        fig.savefig(outdir / f"{prefix}_timeseries_{col}.png", dpi=300)
        plt.close(fig)


def compute_decomp_dd(
    wt_decomp_csv: Path,
    mut_decomp_csv: Path,
) -> pd.DataFrame:
    """
    Compute per-residue ΔΔE = MUT - WT for decomposition CSVs.
    Returns table with residue_id + ΔΔ for all numeric columns.
    """
    wt = load_decomp_perframe_csv(wt_decomp_csv)
    mut = load_decomp_perframe_csv(mut_decomp_csv)


    rid_wt = guess_residue_id_column(wt)
    rid_mut = guess_residue_id_column(mut)

    # Standardize residue id column name
    wt = wt.rename(columns={rid_wt: "residue_id"})
    mut = mut.rename(columns={rid_mut: "residue_id"})

    # Identify numeric cols intersection
    wt_num = [c for c in wt.columns if c != "residue_id" and pd.api.types.is_numeric_dtype(wt[c])]
    mut_num = [c for c in mut.columns if c != "residue_id" and pd.api.types.is_numeric_dtype(mut[c])]
    common = sorted(set(wt_num) & set(mut_num))
    if not common:
        raise ValueError("No common numeric columns found between WT and MUT decomposition CSVs.")

    merged = wt[["residue_id"] + common].merge(mut[["residue_id"] + common], on="residue_id", suffixes=("_wt", "_mut"))
    out = pd.DataFrame({"residue_id": merged["residue_id"]})
    for c in common:
        out[f"dd_{c}"] = merged[f"{c}_mut"] - merged[f"{c}_wt"]
    return out


def plot_decomp_dd_top(
    dd_df: pd.DataFrame,
    outpng: Path,
    metric_col: str,
    top_n: int = 20,
    title: str = "",
):
    """
    Plot top residues by |ΔΔmetric| (largest magnitude change).
    For "binding gets worse" interpretation: positive ΔΔTOTAL often indicates weaker binding (depending on sign convention).
    """
    if metric_col not in dd_df.columns:
        raise ValueError(f"metric_col '{metric_col}' not in dd_df columns: {list(dd_df.columns)}")

    sub = dd_df[["residue_id", metric_col]].dropna().copy()
    sub["abs"] = sub[metric_col].abs()
    sub = sub.sort_values("abs", ascending=False).head(top_n)

    fig, ax = plt.subplots(figsize=(7.6, 4.4))
    ax.barh(sub["residue_id"].astype(str).values[::-1], sub[metric_col].values[::-1])
    ax.set_xlabel("ΔΔ Energy (kcal/mol)  (MUT − WT)")
    ax.set_title(title or f"Top {top_n} residues by |{metric_col}|")
    style_axes(ax)
    plt.tight_layout()
    fig.savefig(outpng, dpi=300)
    plt.close(fig)


# =============================================================================
# Report generation
# =============================================================================
def run_report(
    root: Path,
    outdir: Path,
    isoforms: Optional[List[str]],
    systems: Optional[List[str]],
    wt_regex: str,
    mut_regex: str,
    make_timeseries: bool,
    topn: int,
):
    outdir.mkdir(parents=True, exist_ok=True)

    all_dtotal_rows = []
    all_dd_rows = []

    wt_pat = re.compile(wt_regex)
    mut_pat = re.compile(mut_regex)
    
    for iso_dir in discover_isoforms(root, isoforms):
        iso = iso_dir.name
        sys_dirs = discover_systems(iso_dir, systems)
        if not sys_dirs:
            continue

        # group systems into WT-like and MUT-like by regex
        wt_systems = [p for p in sys_dirs if wt_pat.search(p.name)]
        mut_systems = [p for p in sys_dirs if mut_pat.search(p.name)]

        # per-system processing
        for sys_dir in sys_dirs:
            system = sys_dir.name
            reps = discover_replicates(sys_dir)
            if not reps:
                continue

            # output mirrors input
            base_out = outdir / iso / system
            base_out.mkdir(parents=True, exist_ok=True)

            rep_tables: Dict[str, pd.DataFrame] = {}
            for rep_dir in reps:
                rep = rep_dir.name
                files = expected_files(rep_dir)
                if not files.dat.exists():
                    continue

                # parse components
                ddf = parse_delta_from_dat(files.dat)
                rep_tables[rep] = ddf

                # replicate-level plots
                rep_out = base_out / rep
                rep_out.mkdir(parents=True, exist_ok=True)

                rep_plot_df = ddf.rename(columns={"avg": "mean"})[["component", "mean", "sem"]]
                
                plot_components_bar(
                    rep_plot_df,
                    rep_out / f"{iso}_{system}_{rep}_components.png",
                    title=f"{iso} / {system} / {rep}: Components (mean ± SEM across frames)",
                )
                plot_total_only(
                    rep_plot_df,
                    rep_out / f"{iso}_{system}_{rep}_dtotal.png",
                    title=f"{iso} / {system} / {rep}: ΔTOTAL (mean ± SEM across frames)",
                )

                # optional timeseries
                if make_timeseries and files.csv.exists():
                    plot_timeseries(
                        files.csv,
                        rep_out,
                        prefix=f"{iso}_{system}_{rep}",
                    )

                # collect ΔTOTAL numbers
                dfi = ddf.set_index("component")
                if "ΔTOTAL" in dfi.index:
                    all_dtotal_rows.append({
                        "isoform": iso,
                        "system": system,
                        "replicate": rep,
                        "dtotal_avg": float(dfi.loc["ΔTOTAL", "avg"]),
                        "dtotal_sem_frame": float(dfi.loc["ΔTOTAL", "sem"]) if pd.notna(dfi.loc["ΔTOTAL", "sem"]) else np.nan,
                        "dat": str(files.dat),
                    })

            # system-level summary across replicates
            if len(rep_tables) >= 2:
                sys_sum = summarize_components_across_reps(rep_tables)
                sys_sum.to_csv(base_out / f"{iso}_{system}_summary_across_reps.csv", index=False)

                plot_components_bar(
                    sys_sum,
                    base_out / f"{iso}_{system}_components_mean_across_reps.png",
                    title=f"{iso} / {system}: Components (mean ± SEM across replicates)",
                )
                plot_total_only(
                    sys_sum,
                    base_out / f"{iso}_{system}_dtotal_mean_across_reps.png",
                    title=f"{iso} / {system}: ΔTOTAL (mean ± SEM across replicates)",
                )

        # WT vs MUT comparisons within isoform
        # Strategy: if multiple WT-like or MUT-like systems exist, compare all combinations.
        # Typical case: one WT and one MUT.
        for wt_sys in wt_systems:
            for mut_sys in mut_systems:
                wt_name, mut_name = wt_sys.name, mut_sys.name
                comp_out = outdir / iso / "_comparisons" / f"{mut_name}_minus_{wt_name}"
                comp_out.mkdir(parents=True, exist_ok=True)

                # load replicate component tables for both
                wt_rep_tables = {}
                mut_rep_tables = {}

                for rep_dir in discover_replicates(wt_sys):
                    rep = rep_dir.name
                    files = expected_files(rep_dir)
                    if files.dat.exists():
                        wt_rep_tables[rep] = parse_delta_from_dat(files.dat)

                for rep_dir in discover_replicates(mut_sys):
                    rep = rep_dir.name
                    files = expected_files(rep_dir)
                    if files.dat.exists():
                        mut_rep_tables[rep] = parse_delta_from_dat(files.dat)

                if not wt_rep_tables or not mut_rep_tables:
                    continue

                # paired ΔΔ components
                dd_comp = paired_delta_delta_components(wt_rep_tables, mut_rep_tables)
                dd_comp.to_csv(comp_out / f"{iso}_{mut_name}_minus_{wt_name}_dd_components.csv", index=False)

                plot_components_bar(
                    dd_comp.rename(columns={"mean_dd": "mean", "sem_dd": "sem"})[["component", "mean", "sem"]],
                    comp_out / f"{iso}_{mut_name}_minus_{wt_name}_dd_components.png",
                    title=f"{iso}: ΔΔ Components (MUT − WT)\n{mut_name} − {wt_name}",
                    ylabel="ΔΔ Energy (kcal/mol)",
                )

                plot_total_only(
                    dd_comp.rename(columns={"mean_dd": "mean", "sem_dd": "sem"})[["component", "mean", "sem"]],
                    comp_out / f"{iso}_{mut_name}_minus_{wt_name}_dd_dtotal.png",
                    title=f"{iso}: ΔΔTOTAL (MUT − WT)\n{mut_name} − {wt_name}",
                    value_col="mean",
                    err_col="sem",
                )

                # store for a global summary table
                dd_idx = dd_comp.set_index("component")
                if "ΔTOTAL" in dd_idx.index:
                    all_dd_rows.append({
                        "isoform": iso,
                        "wt_system": wt_name,
                        "mut_system": mut_name,
                        "dd_dtotal_mean": float(dd_idx.loc["ΔTOTAL", "mean_dd"]),
                        "dd_dtotal_sem": float(dd_idx.loc["ΔTOTAL", "sem_dd"]),
                    })

                # paired per-residue ΔΔ decomposition (top N) if decomp CSV exists
                # We compute dd per replicate and then average dd across replicates.
                dd_decomp_tables = []
                common_reps = sorted(set(discover_replicates(wt_sys_i := wt_sys)) and set())
                # build actual common reps
                wt_reps = {p.name: p for p in discover_replicates(wt_sys)}
                mut_reps = {p.name: p for p in discover_replicates(mut_sys)}
                common = sorted(set(wt_reps) & set(mut_reps))

                for rep in common:
                    wt_files = expected_files(wt_reps[rep])
                    mut_files = expected_files(mut_reps[rep])
                    if wt_files.decomp_csv.exists() and mut_files.decomp_csv.exists():
                        dd_df = compute_decomp_dd(wt_files.decomp_csv, mut_files.decomp_csv)
                        dd_df["replicate"] = rep
                        dd_decomp_tables.append(dd_df)

                if dd_decomp_tables:
                    big = pd.concat(dd_decomp_tables, ignore_index=True)

                    # average ΔΔ per residue across replicates (paired)
                    # choose a metric col for ranking
                    # prefer dd_TOTAL (or dd_Delta_TOTAL) if present
                    dd_cols = [c for c in big.columns if c.startswith("dd_")]
                    rid = "residue_id"
                    if dd_cols:
                        # aggregate mean across replicates
                        agg = big.groupby(rid)[dd_cols].mean().reset_index()

                        # pick a total-like dd column for ranking
                        total_like = None
                        for c in dd_cols:
                            if c.lower().endswith("total"):
                                total_like = c
                                break
                        if total_like is None:
                            total_like = dd_cols[0]

                        agg.to_csv(comp_out / f"{iso}_{mut_name}_minus_{wt_name}_dd_decomp_by_residue.csv", index=False)

                        plot_decomp_dd_top(
                            agg,
                            comp_out / f"{iso}_{mut_name}_minus_{wt_name}_dd_decomp_top{topn}.png",
                            metric_col=total_like,
                            top_n=topn,
                            title=f"{iso}: per-residue ΔΔ (MUT − WT)\nTop {topn} by |{total_like}|\n{mut_name} − {wt_name}",
                        )

    # write global summary tables
    if all_dtotal_rows:
        pd.DataFrame(all_dtotal_rows).to_csv(outdir / "all_isoforms_systems_replicates_dtotal.csv", index=False)
    if all_dd_rows:
        pd.DataFrame(all_dd_rows).to_csv(outdir / "all_isoforms_dd_dtotal_summary.csv", index=False)


# =============================================================================
# CLI
# =============================================================================
def main():
    ap = argparse.ArgumentParser(
        description="Generate local publication-style plots from gmx_MMPBSA outputs across isoforms/systems/replicates, including paired WT vs MUT comparisons."
    )
    ap.add_argument("--root", required=True, help="Path to local gmx_MMPBSA root (contains isoform dirs)")
    ap.add_argument("--outdir", default="mmpbsa_figures", help="Output directory")
    ap.add_argument("--isoforms", nargs="*", default=None, help="Subset isoforms, e.g., smb2 sma")
    ap.add_argument("--systems", nargs="*", default=None, help="Subset systems, e.g., wt_atp r108w_atp")
    ap.add_argument("--wt-regex", default=r"^wt_", help="Regex to identify WT systems (default: '^wt_')")
    ap.add_argument("--mut-regex", default=r"^(?!wt_).+", help="Regex to identify MUT systems (default: 'not wt_')")
    ap.add_argument("--timeseries", action="store_true", help="Also plot per-frame time series from -eo CSV")
    ap.add_argument("--topn", type=int, default=20, help="Top N residues for ΔΔ decomposition plots")
    args = ap.parse_args()

    root = Path(args.root).expanduser().resolve()
    outdir = Path(args.outdir).expanduser().resolve()

    run_report(
        root=root,
        outdir=outdir,
        isoforms=args.isoforms,
        systems=args.systems,
        wt_regex=args.wt_regex,
        mut_regex=args.mut_regex,
        make_timeseries=args.timeseries,
        topn=args.topn,
    )

    print(f"Done. Figures and tables written to: {outdir}")


if __name__ == "__main__":
    main()

