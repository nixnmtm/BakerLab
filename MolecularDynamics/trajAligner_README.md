# TrajAligner

**TrajAligner** is a Python command‑line utility for aligning molecular
dynamics trajectories using **MDAnalysis**, with support for:

-   Direct atom selections
-   Paired mobile/reference selections
-   **Sequence‑mapped residue alignment** for homologous proteins with
    different residue numbering

This tool is particularly useful for comparative MD studies such as
**Ref_Protein vs Mob_Protein**, where equivalent structural regions exist but residue
numbering differs.

------------------------------------------------------------------------

# Features

-   Align a **mobile trajectory** to a **reference structure or
    trajectory**
-   Three fitting modes:
    -   Direct selection
    -   Paired mobile/reference selections
    -   Sequence‑mapped alignment
-   Optional **strict mode** for sequence-mapped fitting
-   Optional **dry‑run mode** for debugging selections
-   Exportable:
    -   final selection strings
    -   mapped residue tables

------------------------------------------------------------------------

# Why Sequence‑Mapped Alignment?

In comparative MD studies:

-   Homologous proteins often have **different residue numbering**
-   Insertions/deletions shift residue indices
-   Direct selections like `resid 100-200` may not correspond across
    systems

Example:

  Reference   Mobile
  ----------- --------
  Ref_Protein        Mob_Protein

Residue **120 in Ref_Protein** may correspond to **118 in Mob_Protein**.

TrajAligner automatically determines this mapping using **sequence
alignment**.

------------------------------------------------------------------------

# Requirements

Python ≥ 3.9

Required packages:

    MDAnalysis
    Biopython

Install with:

``` bash
pip install MDAnalysis biopython
```

------------------------------------------------------------------------

# Input Concepts

## Mobile System

System whose trajectory will be aligned.

    --mob-struct
    --mob-traj

## Reference System

    --ref-struct
    --ref-traj (optional)

Usually only the reference structure is required.

------------------------------------------------------------------------

# Selection Modes

TrajAligner supports **three alignment modes**.

Only **one mode should be used per run**.

------------------------------------------------------------------------

# 1. Direct Selection Mode

Same selection applied to both systems.

Example:

``` bash
--select "protein and name CA"
```

Recommended when:

-   systems share similar topology
-   residue numbering differences are irrelevant

------------------------------------------------------------------------

# 2. Paired Selection Mode

Specify different selections for reference and mobile.

Example:

``` bash
--select-paired \
--select-mob "protein and resid 10:200 and name CA" \
--select-ref "protein and resid 15:205 and name CA"
```

Use when you already know exact equivalent residues.

------------------------------------------------------------------------

# 3. Sequence‑Mapped Mode

Automatically map reference residues to mobile residues using sequence
alignment.

Example:

``` bash
--select-seqmap \
--seq-map on \
--ref-resids 1-800 \
--fit-atoms backbone
```

Recommended for homologous proteins with different numbering.

------------------------------------------------------------------------

# Command Line Arguments

## Input / Output

### Mobile structure

    --mob-struct

Examples:

-   `.gro`
-   `.pdb`
-   `.psf`
-   `.tpr`

### Mobile trajectory

    --mob-traj

Examples:

-   `.xtc`
-   `.trr`
-   `.dcd`

### Reference structure

    --ref-struct

If omitted, the mobile structure is used.

### Reference trajectory

    --ref-traj

Optional. Usually unnecessary.

### Reference frame

    --ref-frame

Default:

    0

### Output trajectory

    -o / --out

Default:

    traj.aligned.xtc

------------------------------------------------------------------------

# Fitting Options

### Direct selection

    --select "protein and name CA"

------------------------------------------------------------------------

### Paired selection

    --select-paired
    --select-mob
    --select-ref

------------------------------------------------------------------------

### Sequence mapped mode

    --select-seqmap

------------------------------------------------------------------------

### Sequence mapping behavior

    --seq-map {off,on,auto}

  option   description
  -------- --------------------------------------------
  off      never use mapping
  on       always use mapping
  auto     fallback to mapping if selections mismatch

------------------------------------------------------------------------

### Segment IDs

    --segid-ref
    --segid-mob

Used when structures contain multiple chains.

------------------------------------------------------------------------

### Fit atoms

    --fit-atoms {CA, backbone, custom}

Options:

  option     meaning
  ---------- ------------------------
  CA         C‑alpha atoms
  backbone   backbone atoms
  custom     user defined selection

Custom example:

    --fit-atoms custom --fit-atom-selection "name CA or name CB"

------------------------------------------------------------------------

### Reference residues

    --ref-resids

Example:

    --ref-resids 10-50,60,72-90

Residues are interpreted in **reference numbering** and mapped
automatically.

------------------------------------------------------------------------

### Strict mode

    --strict

Raises an error if:

-   residues cannot be mapped
-   atoms are missing
-   selections mismatch

Without strict mode problematic residues are skipped.

------------------------------------------------------------------------

# Diagnostic Options

### Verbose mode

    -v

### Dry run

    --dry-run

Shows mapping and selections without writing output trajectory.

### Export selections

    --sel-out selection.txt

### Export mapping

    --map-out mapping.tsv

Example mapping file:

    ref_resid    mob_resid
    120          118
    121          119

------------------------------------------------------------------------

# Examples

## Standard CA alignment

``` bash
python trajAligner.py \
  --mob-struct conf.gro \
  --mob-traj traj.xtc \
  --ref-struct ref.gro \
  --select "protein and name CA" \
  -o traj.aligned.xtc
```

------------------------------------------------------------------------

## Paired selections

``` bash
python trajAligner.py \
  --mob-struct Mob_Protein.psf \
  --mob-traj Mob_Protein.xtc \
  --ref-struct Ref_Protein.psf \
  --select-paired \
  --select-mob "protein and name CA" \
  --select-ref "protein and name CA" \
  -o aligned.xtc
```

------------------------------------------------------------------------

## Sequence mapped alignment

``` bash
python trajAligner.py \
  --mob-struct Mob_Protein.psf \
  --mob-traj prod.xtc \
  --ref-struct Ref_Protein.psf \
  --select-seqmap \
  --seq-map on \
  --fit-atoms CA \
  -o Mob_Protein_aligned.xtc
```

------------------------------------------------------------------------

## Residue range mapped alignment

``` bash
python trajAligner.py \
  --mob-struct Mob_Protein.psf \
  --mob-traj prod.xtc \
  --ref-struct Ref_Protein.psf \
  --select-seqmap \
  --seq-map on \
  --ref-resids 1-800 \
  --fit-atoms backbone \
  --dry-run \
  --sel-out selection.out \
  --map-out mapping.out
```

------------------------------------------------------------------------

# Typical Ref_Protein vs Mob_Protein Workflow

1.  Define residues in **Ref_Protein numbering**
2.  Run sequence‑mapped alignment
3.  Export mapping table
4.  Verify selections using `--dry-run`
5.  Perform full trajectory alignment

------------------------------------------------------------------------

# Troubleshooting

### Atom mismatch error

Try:

    --fit-atoms CA

instead of backbone.

------------------------------------------------------------------------

### No matched residues

Possible causes:

-   wrong `segid`
-   incorrect residue list
-   incomplete structures

------------------------------------------------------------------------

### Unexpected mapping

Inspect:

    mapping.out
    selection.out

------------------------------------------------------------------------

# Output Files

  file            description
  --------------- -------------------------
  aligned.xtc     aligned trajectory
  selection.out   final selection strings
  mapping.out     residue mapping table

------------------------------------------------------------------------

# Dependencies

This tool relies on:

-   **MDAnalysis**
-   **Biopython**

Please cite these packages if used in publications.

------------------------------------------------------------------------

# License

GPL v3

------------------------------------------------------------------------
