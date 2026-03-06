#!/usr/bin/env python3
"""
Align trajectories with optional sequence-mapped residue selection.

@author:     Nixon Raj
@copyright:  2023 Nationwide Childrens Hospital, Abigail Wexner Research Institute
@license:    GPL v. 3
@contact:    nixon.raj@nationwidechildrens.org
"""

import argparse
import os
import os.path as path
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Union

import MDAnalysis as mda
import MDAnalysis.analysis.align as mdalign

PROGRAM_VERSION = "v3.1"


# CLI parsing

def parse_args():
    parser = argparse.ArgumentParser(
        description="Align trajectories (supports direct selection, paired selection, or sequence-mapped selection).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-V", "--version", action="version", version=PROGRAM_VERSION)

    beh = parser.add_argument_group("Behavior / Diagnostics")
    beh.add_argument("-v", "--verbose", action="store_true", help="Verbose logging.")
    beh.add_argument("--dry-run", action="store_true",
                     help="Do everything except writing the aligned trajectory. Prints final selections and counts.")
    beh.add_argument("--sel-out", default=None,
                     help="Write final mobile/ref selection strings (and counts) to this file.")
    beh.add_argument("--map-out", default=None,
                     help="Write resid mapping table (ref_resid -> mob_resid) used for seq-map selection to TSV.")

    io = parser.add_argument_group("Input / Output")
    io.add_argument("--mob-struct", required=True, help="Mobile structure file (e.g., .gro/.pdb/.tpr).")
    io.add_argument("--mob-traj", required=True, help="Mobile trajectory file (e.g., .xtc/.trr).")
    io.add_argument("--ref-struct", default=None,
                    help="Reference structure file (PSF, PDB, GRO). If omitted, ref-struct = mob-struct.")
    io.add_argument("--ref-traj", default=None,
                    help=(
                        "Optional reference trajectory (can be a static PDB or GRO).\n"
                        "If PSF is provided as structure we need any coordinate file "
                        "(PDB, CRD, or XTC)."))
    io.add_argument("--ref-frame", type=int, default=0,
                    help="Reference frame index to use (only relevant if --ref-traj is provided).")
    io.add_argument("-o", "--out", dest="trajout",
                    default=path.join(os.curdir, "traj.aligned.xtc"),
                    help="Output aligned trajectory path.")

    sel = parser.add_argument_group("Fitting / Selection (choose ONE mode)")
    mode = sel.add_mutually_exclusive_group()

    # Mode A: one selection string applied to both universes
    mode.add_argument("--select", default=None,
                      help='Direct selection applied to BOTH ref and mobile (e.g., "protein and name CA").')

    # Mode B: explicit paired selections
    mode.add_argument("--select-paired", action="store_true",
                      help="Use paired selections: --select-mob and --select-ref.")

    sel.add_argument("--select-mob", default=None,
                     help='Mobile selection string (paired mode only).')
    sel.add_argument("--select-ref", default=None,
                     help='Reference selection string (paired mode only).')

    # Mode C: sequence-mapped residue selection
    mode.add_argument("--select-seqmap", action="store_true",
                      help="Use sequence-mapped residue selection (ref resid list mapped to mobile via alignment).")
    sel.add_argument("--seq-map", choices=["off", "on", "auto"], default="off",
                     help=("off: never use sequence mapping; on: always use sequence mapping when seqmap mode or when needed; "
                           "auto: try direct/paired selection first, fall back to seq mapping only if paired atom counts mismatch."))

    sel.add_argument("--segid-ref", default=None,
                     help="segid for reference protein segment/chain (optional).")
    sel.add_argument("--segid-mob", default=None,
                     help="segid for mobile protein segment/chain (optional).")

    sel.add_argument("--fit-atoms", choices=["CA", "backbone", "custom"], default="CA",
                     help="Which atoms within residues to use for fitting in seq-map mode.")
    sel.add_argument("--fit-atom-selection", default=None,
                     help='If --fit-atoms custom, provide MDAnalysis selection (e.g., "name CA or name CB").')

    sel.add_argument("--ref-resids", default=None,
                     help="Comma/range list in REF resid numbering (e.g. 10-50,60,72-90). Only used in seq-map mode.")

    sel.add_argument("--strict", action="store_true",
                     help=("Strict mode for seq-map: if a requested ref resid cannot be mapped or lacks atoms in either system, "
                           "raise an error instead of skipping."))

    args = parser.parse_args()

    # Defaults: if user did not select a mode, treat as Mode A with default selection
    if (args.select is None) and (not args.select_paired) and (not args.select_seqmap):
        args.select = "protein and name CA"

    # Validate paired mode inputs
    if args.select_paired:
        if not args.select_mob or not args.select_ref:
            parser.error("--select-paired requires both --select-mob and --select-ref.")

    # Validate seqmap mode inputs
    if args.select_seqmap:
        if args.seq_map == "off":
            # Still allow, but it's contradictory; make it explicit
            parser.error("--select-seqmap requires --seq-map on (or auto). Set --seq-map on.")
        if args.fit_atoms == "custom" and not args.fit_atom_selection:
            parser.error("--fit-atoms custom requires --fit-atom-selection.")

    # Ref defaults
    if args.ref_struct is None:
        args.ref_struct = args.mob_struct

    return args

# Utility: residue list parsing
def parse_resid_list(s: Optional[str]) -> Optional[List[int]]:
    """Parse '10-12,20,30-33' -> sorted unique ints."""
    if s is None:
        return None
    out = set()
    for part in s.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-", 1)
            a, b = int(a), int(b)
            lo, hi = (a, b) if a <= b else (b, a)
            out.update(range(lo, hi + 1))
        else:
            out.add(int(part))
    return sorted(out)

# Sequence mapping utilities
def _validate_segid(u: mda.Universe, segid: Optional[str], label: str):
    if segid:
        if u.select_atoms(f"segid {segid} and protein").n_atoms == 0:
            raise ValueError(f"{label} segid '{segid}' has no protein atoms. Check segid names in your structure.")

def structure_sequence_and_resids(u: mda.Universe, segid: Optional[str] = None) -> Tuple[str, List[int]]:
    """
    Build 1-letter AA sequence and corresponding resid list from the structure.
    Uses protein residues in order.
    """
    from MDAnalysis.lib.util import convert_aa_code

    sel = "protein" if not segid else f"segid {segid} and protein"
    residues = u.select_atoms(sel).residues

    seq = []
    resids = []
    for r in residues:
        try:
            seq.append(convert_aa_code(r.resname))  # 3->1
            resids.append(int(r.resid))
        except Exception:
            # skip non-standard/unknown residues
            continue

    return "".join(seq), resids

def global_sequence_align(seqA: str, seqB: str) -> Tuple[str, str, float]:
    """
    Returns aligned strings (with '-' gaps). Uses Biopython.
    """
    try:
        # Prefer modern aligner if available
        from Bio.Align import PairwiseAligner
        aligner = PairwiseAligner()
        aligner.mode = "global"
        # Scoring roughly similar to original pairwise2 settings
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -0.5

        aln = next(iter(aligner.align(seqA, seqB)))
        # Reconstruct gapped strings
        # PairwiseAligner stores aligned coordinates; we convert to gapped strings
        # This is a compact, reliable conversion:
        aA = []
        aB = []
        iA = iB = 0
        for (a0, a1), (b0, b1) in zip(aln.aligned[0], aln.aligned[1]):
            # gaps before this block
            while iA < a0 and iB < b0:
                aA.append(seqA[iA]); aB.append(seqB[iB]); iA += 1; iB += 1
            while iA < a0:
                aA.append(seqA[iA]); aB.append("-"); iA += 1
            while iB < b0:
                aA.append("-"); aB.append(seqB[iB]); iB += 1
            # aligned block
            while iA < a1 and iB < b1:
                aA.append(seqA[iA]); aB.append(seqB[iB]); iA += 1; iB += 1
        # tail
        while iA < len(seqA) and iB < len(seqB):
            aA.append(seqA[iA]); aB.append(seqB[iB]); iA += 1; iB += 1
        while iA < len(seqA):
            aA.append(seqA[iA]); aB.append("-"); iA += 1
        while iB < len(seqB):
            aA.append("-"); aB.append(seqB[iB]); iB += 1

        return "".join(aA), "".join(aB), float(aln.score)

    except Exception:
        # Fallback to pairwise2 if needed
        try:
            from Bio import pairwise2
        except ImportError as e:
            raise ImportError("Biopython not found. Install with: pip install biopython") from e

        aln = pairwise2.align.globalms(seqA, seqB, 2, -1, -10, -0.5)[0]
        return aln[0], aln[1], float(aln[2])

def build_resid_map_from_alignment(alnA: str, alnB: str, residsA: List[int], residsB: List[int]) -> Dict[int, int]:
    """
    Map residA -> residB for alignment columns where both are non-gap.
    """
    iA = 0
    iB = 0
    m = {}
    for cA, cB in zip(alnA, alnB):
        if cA != "-" and cB != "-":
            m[residsA[iA]] = residsB[iB]
        if cA != "-":
            iA += 1
        if cB != "-":
            iB += 1
    return m

def _fit_atomsel(fit_atoms: str, fit_atom_selection: Optional[str]) -> str:
    if fit_atoms == "CA":
        return "name CA"
    if fit_atoms == "backbone":
        return "backbone"
    if fit_atoms == "custom":
        return fit_atom_selection
    raise ValueError(f"Unknown fit_atoms: {fit_atoms}")

def build_seqmap_selection(
    ref_u: mda.Universe,
    mob_u: mda.Universe,
    segid_ref: Optional[str],
    segid_mob: Optional[str],
    fit_atomsel: str,
    ref_resids_keep: Optional[List[int]],
    strict: bool,
    verbose: bool,
) -> Tuple[Tuple[str, str], Dict[int, int], List[int], List[int], float]:
    """
    Returns:
      select_tuple = (mob_sel_str, ref_sel_str)
      resid_map (full mapping)
      used_ref_resids, used_mob_resids (actually used after filtering)
      alignment_score
    """
    _validate_segid(ref_u, segid_ref, "ref")
    _validate_segid(mob_u, segid_mob, "mob")

    ref_seq, ref_resids = structure_sequence_and_resids(ref_u, segid=segid_ref)
    mob_seq, mob_resids = structure_sequence_and_resids(mob_u, segid=segid_mob)

    if not ref_seq or not mob_seq:
        raise ValueError("Could not derive protein sequence from structure(s). Check protein selection and segid.")

    aln_ref, aln_mob, score = global_sequence_align(ref_seq, mob_seq)
    resid_map = build_resid_map_from_alignment(aln_ref, aln_mob, ref_resids, mob_resids)

    if ref_resids_keep is None:
        ref_resids_keep = sorted(resid_map.keys())

    base_ref = "protein" if not segid_ref else f"segid {segid_ref} and protein"
    base_mob = "protein" if not segid_mob else f"segid {segid_mob} and protein"

    used_ref = []
    used_mob = []

    for r in ref_resids_keep:
        if r not in resid_map:
            if strict:
                raise ValueError(f"[strict] ref resid {r} is not mappable (gap/unaligned).")
            continue
        rm = resid_map[r]

        a_ref = ref_u.select_atoms(f"{base_ref} and ({fit_atomsel}) and resid {r}")
        a_mob = mob_u.select_atoms(f"{base_mob} and ({fit_atomsel}) and resid {rm}")

        if a_ref.n_atoms == 0 or a_mob.n_atoms == 0:
            if strict:
                raise ValueError(f"[strict] missing fit atoms for resid pair ref={r} mob={rm} using '{fit_atomsel}'.")
            continue
        if a_ref.n_atoms != a_mob.n_atoms:
            if strict:
                raise ValueError(f"[strict] atomcount mismatch for resid pair ref={r} mob={rm}: "
                                 f"ref={a_ref.n_atoms} mob={a_mob.n_atoms}")
            continue

        used_ref.append(r)
        used_mob.append(rm)

    if not used_ref:
        raise ValueError("No matched residues found after filtering. "
                         "Try --fit-atoms CA, check segids, or loosen strict mode.")

    ref_sel = f"{base_ref} and ({fit_atomsel}) and resid " + " ".join(map(str, used_ref))
    mob_sel = f"{base_mob} and ({fit_atomsel}) and resid " + " ".join(map(str, used_mob))

    nref = ref_u.select_atoms(ref_sel).n_atoms
    nmob = mob_u.select_atoms(mob_sel).n_atoms
    if nref != nmob:
        raise ValueError(f"Seq-map selection mismatch: ref={nref}, mob={nmob} (this should not happen).")

    if verbose:
        print(f"[seq-map] alignment score: {score}")
        print(f"[seq-map] mapped positions (total): {len(resid_map)}")
        print(f"[seq-map] residues used (after filtering): {len(used_ref)}")
        print(f"[seq-map] atoms used: {nref}")

    return (mob_sel, ref_sel), resid_map, used_ref, used_mob, score

# Selection resolution
@dataclass
class FinalSelection:
    select: Union[str, Tuple[str, str]]
    mode: str  # "direct" | "paired" | "seqmap" | "auto->seqmap"
    details: str


def resolve_selection(args, ref_u: mda.Universe, mob_u: mda.Universe) -> Tuple[FinalSelection, Optional[Dict[int, int]], Optional[List[int]], Optional[List[int]]]:
    """
    Decide final selection for AlignTraj, honoring mutually exclusive modes and seq-map fallback rules.
    Returns:
      FinalSelection
      resid_map, used_ref_resids, used_mob_resids (only for seq-map modes)
    """
    resid_map = None
    used_ref = used_mob = None

    # If user explicitly chose seqmap mode: always do seqmap (args.seq_map must be on/auto validated earlier)
    if args.select_seqmap:
        fit_atomsel = _fit_atomsel(args.fit_atoms, args.fit_atom_selection)
        ref_keep = parse_resid_list(args.ref_resids)

        select_tuple, resid_map, used_ref, used_mob, _score = build_seqmap_selection(
            ref_u, mob_u,
            segid_ref=args.segid_ref,
            segid_mob=args.segid_mob,
            fit_atomsel=fit_atomsel,
            ref_resids_keep=ref_keep,
            strict=args.strict,
            verbose=args.verbose,
            )
        return FinalSelection(select=select_tuple, mode="seqmap", details=f"fit_atomsel={fit_atomsel}"), resid_map, used_ref, used_mob

    # Mode B: paired direct selection
    if args.select_paired:
        select_tuple = (args.select_mob, args.select_ref)

        # If seq_map is auto: check mismatch and fall back to seqmap if needed
        mob_n = mob_u.select_atoms(select_tuple[0]).n_atoms
        ref_n = ref_u.select_atoms(select_tuple[1]).n_atoms

        if mob_n != ref_n:
            if args.seq_map == "auto":
                fit_atomsel = _fit_atomsel(args.fit_atoms, args.fit_atom_selection)
                ref_keep = parse_resid_list(args.ref_resids)

                select_tuple, resid_map, used_ref, used_mob, _score = build_seqmap_selection(
                    ref_u, mob_u,
                    segid_ref=args.segid_ref,
                    segid_mob=args.segid_mob,
                    fit_atomsel=fit_atomsel,
                    ref_resids_keep=ref_keep,
                    strict=args.strict,
                    verbose=args.verbose,
                )
                return FinalSelection(select=select_tuple, mode="auto->seqmap",
                                     details=f"paired mismatch (mob={mob_n}, ref={ref_n}); fallback fit_atomsel={fit_atomsel}"), resid_map, used_ref, used_mob
            raise ValueError(f"Paired selection atom mismatch: mobile={mob_n}, ref={ref_n}. "
                             f"Fix selections or use --seq-map auto/on with --select-seqmap.")

        return FinalSelection(select=select_tuple, mode="paired", details=f"mob_atoms={mob_n} ref_atoms={ref_n}"), None, None, None

    # Mode A: direct selection applied to both
    if args.select is not None:
        # For direct selection, AlignTraj will apply it to both.
        # If user wants mapping, they should use seqmap mode or paired mode with auto fallback.
        return FinalSelection(select=args.select, mode="direct", details=""), None, None, None

    # Should never get here due to default select
    raise RuntimeError("Internal error: no selection mode resolved.")



# Main function
def align_traj(args):
    # Load reference universe efficiently:
    if args.ref_traj:
        ref_u = mda.Universe(args.ref_struct, args.ref_traj)
        if args.ref_frame < 0 or args.ref_frame >= len(ref_u.trajectory):
            raise ValueError(f"--ref-frame {args.ref_frame} out of range for reference trajectory.")
        ref_u.trajectory[args.ref_frame]
    else:
        ref_u = mda.Universe(args.ref_struct)

    # Mobile universe always needs trajectory
    mob_u = mda.Universe(args.mob_struct, args.mob_traj)

    final_sel, resid_map, used_ref, used_mob = resolve_selection(args, ref_u, mob_u)

    # Diagnostics
    if args.verbose or args.dry_run:
        if isinstance(final_sel.select, tuple):
            mob_sel, ref_sel = final_sel.select
            mob_n = mob_u.select_atoms(mob_sel).n_atoms
            ref_n = ref_u.select_atoms(ref_sel).n_atoms
            print(f"[selection] mode={final_sel.mode} ({final_sel.details})")
            print(f"[selection] mob_n={mob_n} ref_n={ref_n}")
        else:
            # direct selection
            # not guaranteed to match counts (and doesn't need to) because selection is applied separately
            print(f"[selection] mode={final_sel.mode} select='{final_sel.select}'")

    # Write selection info if requested
    if args.sel_out:
        with open(args.sel_out, "w") as fh:
            fh.write(f"mode\t{final_sel.mode}\n")
            fh.write(f"details\t{final_sel.details}\n")
            if isinstance(final_sel.select, tuple):
                mob_sel, ref_sel = final_sel.select
                fh.write(f"mob_select\t{mob_sel}\n")
                fh.write(f"ref_select\t{ref_sel}\n")
                fh.write(f"mob_atoms\t{mob_u.select_atoms(mob_sel).n_atoms}\n")
                fh.write(f"ref_atoms\t{ref_u.select_atoms(ref_sel).n_atoms}\n")
            else:
                fh.write(f"select\t{final_sel.select}\n")

    # Write mapping if requested (only meaningful for seqmap modes)
    if args.map_out and resid_map and used_ref and used_mob:
        with open(args.map_out, "w") as fh:
            fh.write("ref_resid\tmob_resid\n")
            for r, m in zip(used_ref, used_mob):
                fh.write(f"{r}\t{m}\n")

    if args.dry_run:
        print("[dry-run] Not writing output trajectory.")
        return
    # Run alignment
    aligned = mdalign.AlignTraj(
        mob_u,
        ref_u,
        select=final_sel.select,
        filename=args.trajout,
        verbose=args.verbose
    )
    aligned.run()


def main():
    args = parse_args()
    align_traj(args)


if __name__ == "__main__":
    main()                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    