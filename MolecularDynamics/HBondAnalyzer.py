# hbond_analyzer.py
from __future__ import annotations
import os
from collections import defaultdict
from typing import Dict, List, Tuple, Optional, Any

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
from concurrent.futures import ProcessPoolExecutor, as_completed

# ------------------------- low-level helpers -------------------------
def _per_frame_sets_from_results(u: mda.Universe, H: HBA, analyzed_frames: List[int]):
    """Return [set((donor_idx, acceptor_idx)), ...] in analyzed frame order."""
    rows = getattr(getattr(H, "results", None), "hbonds", None)  # ndarray: frame, donor, hydrogen, acceptor, dist, ang
    per_frame = [set() for _ in analyzed_frames]
    if rows is not None:
        m = {}
        for fr, d, h, a, *_ in rows:
            m.setdefault(int(fr), set()).add((int(d), int(a)))
        for i, fi in enumerate(analyzed_frames):
            per_frame[i] = m.get(int(fi), set())
        return per_frame

    # fallback to timeseries styles
    ts = H.timeseries
    if isinstance(ts, list) and len(ts) == len(analyzed_frames):
        out = []
        for fr in ts:
            s = set()
            for row in fr:
                if len(row) >= 6:
                    _, d, _, a, *_ = row
                elif len(row) == 4:
                    d, a, *_ = row
                else:
                    continue
                s.add((int(d), int(a)))
            out.append(s)
        return out

    by_frame = {}
    for row in ts:
        if len(row) >= 6:
            fi, d, _, a, *_ = row
            by_frame.setdefault(int(fi), set()).add((int(d), int(a)))
    return [by_frame.get(int(fi), set()) for fi in analyzed_frames]


def _runlengths_with_censor(b: np.ndarray):
    """Return (lengths_in_frames, event_observed) for each continuous True segment."""
    lengths, observed = [], []
    i, n = 0, b.size
    while i < n:
        if not b[i]:
            i += 1
            continue
        j = i
        while j < n and b[j]:
            j += 1
        lengths.append(j - i)
        observed.append(j < n)  # False if we hit array end → censored
        i = j
    return np.asarray(lengths, int), np.asarray(observed, bool)


def _km_with_censor(times, events):
    """Kaplan–Meier estimator with right-censoring. `times` in same units."""
    if len(times) == 0:
        return np.array([0.0]), np.array([1.0])
    order = np.argsort(times)
    t = np.asarray(times)[order]
    e = np.asarray(events, bool)[order]
    uniq = np.unique(t)
    n_risk = len(t)
    T = [0.0]; S = [1.0]
    for ti in uniq:
        at_t = (t == ti)
        d_i = int(e[at_t].sum())           # events
        c_i = int((~e[at_t]).sum())        # censors
        if n_risk > 0 and d_i > 0:
            S.append(S[-1] * (1.0 - d_i / n_risk))
        else:
            S.append(S[-1])
        T.append(float(ti))
        n_risk -= (d_i + c_i)
    return np.array(T), np.array(S)


def _build_strict_site_selectors(
    u: mda.Universe,
    atp_sel: str = "resname ATP",
    site_scope: str = "protein and around 6 resname ATP",
    include_his_as_acceptor: bool = True,
    include_cys_donor: bool = False,         # set True only if you know Cys is protonated
) -> Tuple[str, str, str, str]:
    """
    Returns (prot_donors, prot_acceptors, atp_donors, atp_acceptors), each intersected with site_scope.

    - Protein donors (heavy atoms only): backbone N; ARG NE/NH1/NH2; LYS NZ; HIS ND1/NE2; TRP NE1;
      SER OG; THR OG1; TYR OH; ASN ND2; GLN NE2; (optional CYS SG).
    - Protein acceptors: backbone O/OXT; ASP OD1/OD2; GLU OE1/OE2; ASN OD1; GLN OE1;
      SER OG; THR OG1; TYR OH; (optional HIS ND1/NE2 as acceptors).
    - ATP donors: N6, O2', O3', O5' (only heavy atoms; HBA finds H’s).
    - ATP acceptors: all O* + adenine N1/N7 (+N3 liberal) + O4'.

    Intersections with `site_scope` keep the search local and fast.
    """
    # ---- protein donors (heavy atoms that can carry H) ----
    donors_core = " or ".join([
        "(backbone and name N)",
        "(resname ARG and (name NE or name NH1 or name NH2))",
        "(resname LYS and name NZ)",
        "(resname HIS and (name ND1 or name NE2))",
        "(resname TRP and name NE1)",
        "(resname SER and name OG)",
        "(resname THR and name OG1)",
        "(resname TYR and name OH)",
        "(resname ASN and name ND2)",
        "(resname GLN and name NE2)",
        "(resname CYS and name SG)" if include_cys_donor else "name NONE"
    ])
    prot_donors = f"({site_scope}) and ({donors_core})"

    # ---- protein acceptors (heavy atoms with lone pairs) ----
    acceptors_list = [
        "(backbone and (name O or name OXT))",
        "(resname ASP and (name OD1 or name OD2))",
        "(resname GLU and (name OE1 or name OE2))",
        "(resname ASN and name OD1)",
        "(resname GLN and name OE1)",
        "(resname SER and name OG)",
        "(resname THR and name OG1)",
        "(resname TYR and name OH)",
    ]
    if include_his_as_acceptor:
        acceptors_list.append("(resname HIS and (name ND1 or name NE2))")
    acceptors_core = " or ".join(acceptors_list)
    prot_acceptors = f"({site_scope}) and ({acceptors_core})"

    # ---- ATP side (use your strict definitions) ----
    atp_donors     = f"({atp_sel}) and (name N6 or name O2' or name O3' or name O5')"
    atp_acceptors  = f"({atp_sel}) and (name O* or name O4' or name N1 or name N7 or name N3)"

    return prot_donors, prot_acceptors, atp_donors, atp_acceptors

# bonds_to_vmd_pairs.py
def _parse_bond_tag(tag: str):
    # "LYS184:NZ -> ATP1:O3G" -> ("LYS",184,"NZ","ATP",1,"O3G")
    left, right = [s.strip() for s in tag.split("->")]
    def split_side(side):
        res, atom = [p.strip() for p in side.split(":")]
        i = 0
        while i < len(res) and res[i].isalpha():
            i += 1
        resn = res[:i]
        resid = int(res[i:])
        return resn, resid, atom
    d_resn, d_resi, d_atom = split_side(left)
    a_resn, a_resi, a_atom = split_side(right)
    return d_resn, d_resi, d_atom, a_resn, a_resi, a_atom

def emit_vmd_pairs(res: dict, top_k: int = 20, min_occ: float = 0.05):
    """
    Print a TCL list HBONDS where each item is:
      {d_resn d_resid d_atom  a_resn a_resid a_atom  occupancy}
    Copy-paste this block into the VMD script below.
    """
    items = sorted(res.items(), key=lambda kv: kv[1].get("occupancy", 0.0), reverse=True)
    rows = []
    for tag, data in items:
        occ = float(data.get("occupancy", 0.0))
        if occ < min_occ:
            continue
        rows.append((*_parse_bond_tag(tag), occ))
        if len(rows) >= top_k:
            break

    print("set HBONDS {")
    for d_resn, d_resi, d_atom, a_resn, a_resi, a_atom, occ in rows:
        # primes in atom names are fine in TCL when inside braces
        print(f"  {{ {d_resn} {d_resi} {{{d_atom}}}  {a_resn} {a_resi} {{{a_atom}}}  {occ:.3f} }}")
    print("}")


# ------------------------- main class -------------------------
class HBondAnalyzer:
    """
    ATP–site H-bond analysis across multiple MD replicas & conditions.
    Core outputs per bond:
      - pooled occupancy across runs
      - per-run occupancy list
      - lifetimes (frames) with right-censoring → Kaplan–Meier S(t)
    """

    # ---- construction ----
    def __init__(
        self,
        step: int = 10,            # analyze every N frames
        d_a_cutoff: float = 3.5,   # Å donor–acceptor distance
        d_h_a_angle: float = 135., # deg D–H–A (looser for base contacts)
        d_h_cutoff: float = 1.2,   # Å to detect D–H pairs
        pad_to_full_window: bool = True
    ):
        self.step = step
        self.d_a = d_a_cutoff
        self.dha = d_h_a_angle
        self.dh = d_h_cutoff
        self.pad = pad_to_full_window

    # ---- selections (auto, strict protein) ----
    def build_protein_ligand_hbond_selections(
        self,
        u: mda.Universe,
        atp_sel: str = "resname ATP",
        scope_sel: str = "protein and around 6 resname ATP",
        include_base_N3: bool = True,
        strict_protein: bool = True,
        include_his_as_acceptor: bool = True,
        include_cys_donor: bool = False,   # set True only if you know CYS is protonated
    ) -> Tuple[str, str, str, str]:
        """
        Build two directional H-bond selections:
          A: protein donors  -> ATP acceptors
          B: ATP donors      -> protein acceptors
    
        If strict_protein=True, constrain protein donors/acceptors to chemically sensible
        heavy atoms (and intersect with `scope_sel`) to avoid spurious pairs like CA, CD, etc.
        """
        atp = u.select_atoms(atp_sel)
    
        # --- helper: only include ATP atoms that actually exist in this topology
        def present(names):
            have = set(at.name for at in atp.atoms)
            return [nm for nm in names if nm in have]
    
        # ---------------- ATP side ----------------
        # Acceptors: all O's + base N1/N7/(optional N3) + O4'
        o_acceptors = [at.name for at in atp if at.name.startswith("O")]
        base = ["N1", "N7"] + (["N3"] if include_base_N3 else [])
        base = present(base)
        atp_acceptor_names = sorted(set(o_acceptors + base + present(["O4'"])))
    
        # Donors: heavy atoms that have hydrogens present (checked in topology)
        donors_candidates = {
            "N6": ["H61", "H62"],
            "O2'": ["H2'", "H2''"],
            "O3'": ["H3'"],
            "O5'": ["H5'", "H5''"],
        }
        have_names = set(at.name for at in atp.atoms)
        atp_donor_names = [
            heavy for heavy, hs in donors_candidates.items()
            if heavy in have_names and any(h in have_names for h in hs)
        ]
    
        def clause(names):
            if not names:
                return "name NONE"
            return " or ".join(f"name {nm}" for nm in names)
    
        atp_acceptors = f"({atp_sel}) and ({clause(atp_acceptor_names)})"
        atp_donors    = f"({atp_sel}) and ({clause(atp_donor_names)})"
    
        # ---------------- protein side ----------------
        if strict_protein:
            # Donors: backbone N; ARG NE/NH1/NH2; LYS NZ; HIS ND1/NE2; TRP NE1;
            #         SER OG; THR OG1; TYR OH; ASN ND2; GLN NE2; (optional CYS SG)
            donors_core = " or ".join([
                "(backbone and name N)",
                "(resname ARG and (name NE or name NH1 or name NH2))",
                "(resname LYS and name NZ)",
                "(resname HIS and (name ND1 or name NE2))",
                "(resname TRP and name NE1)",
                "(resname SER and name OG)",
                "(resname THR and name OG1)",
                "(resname TYR and name OH)",
                "(resname ASN and name ND2)",
                "(resname GLN and name NE2)",
                "(resname CYS and name SG)" if include_cys_donor else "name NONE",
            ])
            prot_donors = f"({scope_sel}) and ({donors_core})"
    
            # Acceptors: backbone O/OXT; ASP OD1/OD2; GLU OE1/OE2; ASN OD1; GLN OE1;
            #           SER OG; THR OG1; TYR OH; (optional HIS ND1/NE2)
            acc_parts = [
                "(backbone and (name O or name OXT))",
                "(resname ASP and (name OD1 or name OD2))",
                "(resname GLU and (name OE1 or name OE2))",
                "(resname ASN and name OD1)",
                "(resname GLN and name OE1)",
                "(resname SER and name OG)",
                "(resname THR and name OG1)",
                "(resname TYR and name OH)",
            ]
            if include_his_as_acceptor:
                acc_parts.append("(resname HIS and (name ND1 or name NE2))")
            prot_acceptors = f"({scope_sel}) and ({' or '.join(acc_parts)})"
        else:
            # original broad behavior (not recommended)
            prot_donors    = f"({scope_sel})"
            prot_acceptors = f"({scope_sel})"
    
        # Direction A: protein donors -> ATP acceptors
        donA, accA = prot_donors, atp_acceptors
        # Direction B: ATP donors    -> protein acceptors
        donB, accB = atp_donors,  prot_acceptors
        return donA, accA, donB, accB

    @staticmethod
    def _hba_one_run(args):
        """Worker: run one directional HBA on a single (top,traj)."""
        (top, traj, start, step, d_a, d_h_a, d_h, donors_sel, acceptors_sel, pad) = args
        import MDAnalysis as mda
        from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
        import numpy as np
        from collections import defaultdict
    
        u = mda.Universe(top, traj)
        dt_ps = float(getattr(u.trajectory, "dt", 1.0))
        dt_eff = dt_ps * step
        analyzed = list(range(start, len(u.trajectory), step))
    
        H = HBA(universe=u,
                donors_sel=donors_sel,
                acceptors_sel=acceptors_sel,
                d_h_cutoff=d_h,
                d_a_cutoff=d_a,
                d_h_a_angle_cutoff=d_h_a,
                update_selections=False)
        H.run(start=start, step=step)
    
        # reuse your helpers:
        per_frame_sets = _per_frame_sets_from_results(u, H, analyzed_frames=analyzed)
        out = defaultdict(lambda: {"present":0, "total":0, "runs":[(dt_eff, len(per_frame_sets))], "bools":[]})
    
        if per_frame_sets:
            all_pairs = set().union(*per_frame_sets)
            idx = {p:i for i,p in enumerate(sorted(all_pairs))}
            B = np.zeros((len(idx), len(per_frame_sets)), dtype=bool)
            for t_idx, present in enumerate(per_frame_sets):
                for p in present:
                    B[idx[p], t_idx] = True
            for p, row in zip(sorted(all_pairs), B):
                d_idx, a_idx = p
                d = u.atoms[d_idx]; a = u.atoms[a_idx]
                tag = f"{d.resname}{int(d.resid)}:{d.name} -> {a.resname}{int(a.resid)}:{a.name}"
                out[tag]["present"] += int(row.sum())
                out[tag]["total"]   += row.size
                out[tag]["bools"].append(row)
    
        return dict(out)

    def _hbond_parallel_over_runs(self, runs, donors_sel, acceptors_sel, start_frames=None, max_workers=3):
        if start_frames is None:
            start_frames = [0]*len(runs)
        tasks = []
        for (top,traj), start in zip(runs, start_frames):
            tasks.append((top, traj, start, self.step, self.d_a, self.dha, self.dh, donors_sel, acceptors_sel, self.pad))
    
        merged = {}
        with ProcessPoolExecutor(max_workers=max_workers) as ex:
            for fut in as_completed([ex.submit(self._hba_one_run, t) for t in tasks]):
                partial = fut.result()
                for tag, d in partial.items():
                    entry = merged.setdefault(tag, {"present":0, "total":0, "bools":[], "dt_eff":None, "tot_frames":[]})
                    entry["present"] += d["present"]
                    entry["total"]   += d["total"]
                    entry["bools"]   += d["bools"]
                    # dt_eff and total frames (assume identical across runs; otherwise store per-run)
                    for dt_eff, nfr in d["runs"]:
                        entry["dt_eff"] = dt_eff
                        entry["tot_frames"].append(nfr)
    
        # finalize to your standard result schema
        results = {}
        for tag, d in merged.items():
            pooled = d["present"]/d["total"] if d["total"] else 0.0
            lengths, events = [], []
            for row in d["bools"]:
                l, e = _runlengths_with_censor(row)
                if l.size:
                    lengths.append(l); events.append(e)
            lengths = np.concatenate(lengths) if lengths else np.array([], int)
            events  = np.concatenate(events)  if events  else np.array([], bool)
            t_ps, S = _km_with_censor(lengths * (d["dt_eff"] or 1.0), events)
            if self.pad and d["tot_frames"]:
                total_ps = max(d["tot_frames"]) * (d["dt_eff"] or 1.0)
                if t_ps[-1] < total_ps:
                    t_ps = np.concatenate([t_ps, [total_ps]])
                    S    = np.concatenate([S,  [S[-1]]])
            results[tag] = {
                "occupancy": pooled,
                "per_run_occupancy": [],  # could compute if you also keep per-run counts
                "lifetimes_frames": lengths,
                "event_observed": events,
                "t": t_ps, "S": S
            }
        return results

    # ---- core analysis over replicas ----
    def _hbond_survival_across_runs(
        self,
        runs: List[Tuple[str, str]],
        donors_sel: Optional[str],
        acceptors_sel: Optional[str],
        hydrogens_sel: Optional[str] = None,
        start_frames: Optional[List[int]] = None,
    ) -> Dict[str, Dict[str, Any]]:
        """Internal engine: one directional pass over a condition's replicas."""
        if start_frames is None:
            start_frames = [0] * len(runs)

        acc = defaultdict(lambda: {
            "present": 0, "total": 0,
            "per_run_occupancy": [],
            "lifetimes_frames": [],
            "event_observed": [],
            "dt_ps_eff": None,
            "total_frames_analyzed": []
        })

        for (top, traj), start in zip(runs, start_frames):
            u = mda.Universe(top, traj)
            dt_ps = float(getattr(u.trajectory, "dt", 1.0))
            dt_eff = dt_ps * self.step
            analyzed = list(range(start, len(u.trajectory), self.step))

            H = HBA(
                universe=u,
                donors_sel=donors_sel,
                hydrogens_sel=hydrogens_sel,
                acceptors_sel=acceptors_sel,
                d_h_cutoff=self.dh,
                d_a_cutoff=self.d_a,
                d_h_a_angle_cutoff=self.dha,
                update_selections=False    # change this to True if 6Ang cutoff has to be dynamic
            )
            H.run(start=start, step=self.step)

            per_frame_sets = _per_frame_sets_from_results(u, H, analyzed_frames=analyzed)
            if not per_frame_sets:
                continue

            all_pairs = set().union(*per_frame_sets)
            if not all_pairs:
                continue

            pair_list = sorted(all_pairs)
            idx_map = {p: i for i, p in enumerate(pair_list)}
            B = np.zeros((len(pair_list), len(per_frame_sets)), dtype=bool)
            for t_idx, present in enumerate(per_frame_sets):
                for p in present:
                    B[idx_map[p], t_idx] = True

            for p, row in zip(pair_list, B):
                d_idx, a_idx = p
                d = u.atoms[d_idx]; a = u.atoms[a_idx]
                tag = f"{d.resname}{int(d.resid)}:{d.name} -> {a.resname}{int(a.resid)}:{a.name}"

                pres = int(row.sum()); tot = row.size
                acc[tag]["present"] += pres
                acc[tag]["total"]   += tot
                acc[tag]["per_run_occupancy"].append(pres / max(tot, 1))

                lengths_frames, observed = _runlengths_with_censor(row)
                if lengths_frames.size:
                    acc[tag]["lifetimes_frames"].append(lengths_frames)
                    acc[tag]["event_observed"].append(observed)

                acc[tag]["dt_ps_eff"] = dt_eff
                acc[tag]["total_frames_analyzed"].append(tot)

        # finalize
        results = {}
        for tag, d in acc.items():
            pooled_occ = (d["present"] / d["total"]) if d["total"] else 0.0
            lengths = np.concatenate(d["lifetimes_frames"]) if d["lifetimes_frames"] else np.array([], int)
            events  = np.concatenate(d["event_observed"])   if d["event_observed"]   else np.array([], bool)
            dt_eff  = d["dt_ps_eff"] if d["dt_ps_eff"] is not None else 1.0
            times_ps = lengths * dt_eff

            t_ps, S = _km_with_censor(times_ps, events)

            if self.pad and d["total_frames_analyzed"]:
                total_ps = max(d["total_frames_analyzed"]) * dt_eff
                if t_ps[-1] < total_ps:
                    t_ps = np.concatenate([t_ps, [total_ps]])
                    S    = np.concatenate([S,  [S[-1]]])

            results[tag] = {
                "occupancy": pooled_occ,
                "per_run_occupancy": d["per_run_occupancy"],
                "lifetimes_frames": lengths,
                "event_observed": events,
                "t": t_ps,   # ps
                "S": S
            }
        return results

    # ---- public: analyze one condition (both directions, merged) ----
    def analyze_condition(
        self,
        runs: List[Tuple[str, str]],
        atp_sel: str = "resname ATP",
        scope_sel: str = "protein and around 6 resname ATP",
        include_base_N3: bool = True,
        start_frames: Optional[List[int]] = None,
        merge_policy: str = "max",  # "max" keep higher-occupancy duplicate; or "sum"
    ) -> Dict[str, Dict[str, Any]]:
        """
        Build selections from topology, run two directional passes, and merge.
        Returns dict[tag] -> metrics.
        """
        u0 = mda.Universe(*runs[0])
        donA, accA, donB, accB = self.build_protein_ligand_hbond_selections(
            u0, atp_sel=atp_sel, scope_sel=scope_sel, include_base_N3=include_base_N3
        )
        res_A = self._hbond_parallel_over_runs(runs, donors_sel=donA, acceptors_sel=accA, start_frames=start_frames)
        res_B = self._hbond_parallel_over_runs(runs, donors_sel=donB, acceptors_sel=accB, start_frames=start_frames)

        merged = dict(res_A)
        for k, v in res_B.items():
            if k not in merged:
                merged[k] = v
            else:
                if merge_policy == "max":
                    if v["occupancy"] > merged[k]["occupancy"]:
                        merged[k] = v
                elif merge_policy == "sum":
                    # combine simple counts/arrays; recompute pooled occupancy & KM approximately
                    merged[k]["per_run_occupancy"] += v["per_run_occupancy"]
                    merged[k]["lifetimes_frames"]   = np.concatenate([merged[k]["lifetimes_frames"], v["lifetimes_frames"]])
                    merged[k]["event_observed"]     = np.concatenate([merged[k]["event_observed"], v["event_observed"]])
                    # re-fit KM
                    t_ps, S = _km_with_censor(merged[k]["t"], merged[k]["event_observed"])
                    merged[k]["t"], merged[k]["S"] = t_ps, S
                    # occupancy: if original totals are unknown, fall back to averaging per-run
                    merged[k]["occupancy"] = np.mean(merged[k]["per_run_occupancy"]) if merged[k]["per_run_occupancy"] else 0.0
        return merged

    # ---- comparison utilities ----
    @staticmethod
    def table_occupancy(res_by_cond: Dict[str, dict], as_percent: bool = True) -> pd.DataFrame:
        """Rows=bonds, Cols=conditions, Values=pooled occupancy."""
        all_bonds = sorted(set().union(*[set(r.keys()) for r in res_by_cond.values()]))
        data = {}
        for cond, res in res_by_cond.items():
            col = {b: res[b]["occupancy"] for b in res}
            data[cond] = [col.get(b, 0.0) for b in all_bonds]
        df = pd.DataFrame(data, index=all_bonds)
        return (100.0 * df) if as_percent else df

    @staticmethod
    def bonds_union(res_by_cond: Dict[str, dict]) -> List[str]:
        return sorted(set().union(*[set(res.keys()) for res in res_by_cond.values()]))

    @staticmethod
    def table_median_survival(res_by_cond: Dict[str, dict], to_ns: bool = True) -> pd.DataFrame:
        """Rows=bonds, Cols=conditions, Values=KM median survival (0 if bond absent)."""
        def median_ps(entry: dict) -> float:
            t, S = entry.get("t", []), entry.get("S", [])
            if t is None or S is None or len(S) == 0:
                return 0.0
            idx = np.where(np.asarray(S) <= 0.5)[0]
            return float(t[idx[0]]) if len(idx) else float(t[-1])
    
        all_bonds = sorted(set().union(*[set(r.keys()) for r in res_by_cond.values()]))
        data = {}
        for cond, res in res_by_cond.items():
            vals = []
            for b in all_bonds:
                med_ps = median_ps(res.get(b, {}))
                vals.append(med_ps / 1000.0 if to_ns else med_ps)
            data[cond] = vals
        units = "ns" if to_ns else "ps"
        df = pd.DataFrame(data, index=all_bonds)
        df.index.name = f"bond_tag"
        df.columns = [f"{c} (median_{units})" for c in df.columns]
        return df

    @staticmethod
    def table_per_run_occupancy_stats(self, res_by_cond: Dict[str, dict], as_percent: bool = True) -> pd.DataFrame:
        """
        Rows=bonds; MultiIndex columns: (condition, stat) where stat ∈ {mean, std, n}.
        Stats are over the list 'per_run_occupancy' within each condition’s res entry.
        """
        bonds = sorted(set().union(*[set(r.keys()) for r in res_by_cond.values()]))
        cols = []
        blocks = []
        for cond, res in res_by_cond.items():
            mean_vals, std_vals, n_vals = [], [], []
            for b in bonds:
                lst = res.get(b, {}).get("per_run_occupancy", [])
                if lst:
                    arr = np.asarray(lst, float)
                    mean_vals.append(arr.mean())
                    std_vals.append(arr.std(ddof=1) if arr.size > 1 else 0.0)
                    n_vals.append(arr.size)
                else:
                    mean_vals.append(0.0); std_vals.append(0.0); n_vals.append(0)
            block = pd.DataFrame({
                (cond, "mean"): mean_vals,
                (cond, "std"):  std_vals,
                (cond, "n"):    n_vals,
            }, index=bonds)
            blocks.append(block)
        df = pd.concat(blocks, axis=1)
        df.index.name = "bond_tag"
        if as_percent:
            # Only scale the mean/std, not n
            df = df.copy()
            for cond in res_by_cond:
                if (cond, "mean") in df.columns:
                    df[(cond, "mean")] = 100.0 * df[(cond, "mean")]
                if (cond, "std") in df.columns:
                    df[(cond, "std")]  = 100.0 * df[(cond, "std")]
        df.columns = pd.MultiIndex.from_tuples(df.columns, names=["condition", "stat"])
        return df

    @staticmethod
    def export_csv_tables(
        self,
        res_by_cond: Dict[str, dict],
        out_dir: str,
        prefix: str = "hbonds",
        include_percent: bool = True,
        include_ns: bool = True
    ) -> Dict[str, str]:
        """
        Write three CSVs to `out_dir`:
          1) occupancy (% or fraction)
          2) median survival (ns or ps)
          3) per-run occupancy stats (mean/std/n; mean/std in % if include_percent)
        Returns a dict of file paths.
        """
        os.makedirs(out_dir, exist_ok=True)
    
        occ = self.table_occupancy(res_by_cond, as_percent=include_percent)
        med = self.table_median_survival(res_by_cond, to_ns=include_ns)
        prs = self.table_per_run_occupancy_stats(res_by_cond, as_percent=include_percent)
    
        units_occ = "pct" if include_percent else "frac"
        units_med = "ns" if include_ns else "ps"
    
        f_occ = os.path.join(out_dir, f"{prefix}_occupancy_{units_occ}.csv")
        f_med = os.path.join(out_dir, f"{prefix}_median_survival_{units_med}.csv")
        f_prs = os.path.join(out_dir, f"{prefix}_perrun_occupancy_stats_{units_occ}.csv")
    
        occ.to_csv(f_occ)
        med.to_csv(f_med)
        prs.to_csv(f_prs)
    
        return {"occupancy": f_occ, "median_survival": f_med, "per_run_stats": f_prs}
    

    
    def top_changed_bonds_by_occupancy(
        self, res_by_cond: Dict[str, dict], k: int = 10, pair: Optional[Tuple[str,str]] = None
    ) -> List[str]:
        df = self.table_occupancy(res_by_cond, as_percent=False)
        if pair is not None:
            a, b = pair
            a_vals = df[a] if a in df.columns else 0.0
            b_vals = df[b] if b in df.columns else 0.0
            delta = (a_vals - b_vals).abs()
        else:
            delta = df.max(axis=1) - df.min(axis=1)
        return list(delta.sort_values(ascending=False).head(k).index)

    def top_changed_bonds_by_median(
        self, res_by_cond: Dict[str, dict], k: int = 10, pair: Optional[Tuple[str,str]] = None
    ) -> List[str]:
        union = self.bonds_union(res_by_cond)
        conds = list(res_by_cond.keys()) if pair is None else list(pair)
        rows = []
        for tag in union:
            if pair is not None:
                a, b = pair
                m_a = self.table_median_survival(res_by_cond[a].get(tag, {}))
                m_b = self.table_median_survival(res_by_cond[b].get(tag, {}))
                delta = abs(m_a - m_b)
            else:
                meds = [self.table_median_survival(res_by_cond[c].get(tag, {})) for c in conds]
                delta = max(meds) - min(meds)
            rows.append((tag, delta))
        rows.sort(key=lambda x: x[1], reverse=True)
        return [tag for tag, _ in rows[:k]]

    # ---- plotting ----
    def _overlay_one_bond_on_axis(
        self, ax, res_by_condition: Dict[str, dict], bond_tag: str,
        to_ns: bool, linewidth: float, alpha: float,
        colors: Optional[Dict[str,str]], cond_order: Optional[List[str]] = None,
        label_with_cond_only: bool = True
    ):
        if cond_order is None:
            cond_order = list(res_by_condition.keys())
        plotted = False
        for cond in cond_order:
            res = res_by_condition[cond]
            if bond_tag not in res:
                continue
            t, S = res[bond_tag]["t"], res[bond_tag]["S"]
            x = (t/1000.0) if to_ns else t
            clr = (colors or {}).get(cond, None)
            lbl = cond if label_with_cond_only else f"{bond_tag} [{cond}]"
            ax.step(x, S, where="post", linewidth=linewidth, alpha=alpha, color=clr, label=lbl)
            plotted = True
        return plotted

    def plot_survival_overlay(
        self, res_by_condition: Dict[str, dict], bond_tag: str, to_ns: bool = True,
        linewidth: float = 2.5, alpha: float = 0.9, figsize=(7,5),
        title: Optional[str] = None, colors: Optional[Dict[str, str]] = None,
        legend_loc: str = "best", ax=None, cond_order: Optional[List[str]] = None,
        label_with_cond_only: bool = True
    ):
        created = False
        if ax is None:
            plt.figure(figsize=figsize); ax = plt.gca(); created = True
        if not self._overlay_one_bond_on_axis(ax, res_by_condition, bond_tag, to_ns, linewidth, alpha, colors, cond_order, label_with_cond_only):
            raise ValueError(f"Bond '{bond_tag}' not found in any provided condition.")
        ax.set_xlabel("Lifetime (ns)" if to_ns else "Lifetime (ps)")
        ax.set_ylabel("Survival S(t)")
        ax.set_ylim(0, 1.05); ax.margins(x=0); ax.grid(True, alpha=0.2)
        ax.legend(frameon=False, loc=legend_loc)
        ax.set_title(title or f"HBond survival — {bond_tag}")
        if created:
            plt.tight_layout(); plt.show()

    def panel_survival_overlays(
        self, res_by_cond: Dict[str, dict], bond_tags: List[str],
        rows: int = 2, cols: int = 3, bonds_per_panel: int = 3,
        colors: Optional[Dict[str,str]] = None, to_ns: bool = True,
        linewidth: float = 2.0, alpha: float = 0.95, sharey: bool = True,
        suptitle: Optional[str] = None, legend_loc: str = "lower left",
        cond_order: Optional[List[str]] = None, label_with_cond_only: bool = True
    ):
        def chunked(seq, n):
            for i in range(0, len(seq), n):
                yield seq[i:i+n]

        n_panels = rows * cols
        tags = bond_tags[:n_panels * bonds_per_panel]
        fig, axes = plt.subplots(rows, cols, figsize=(cols*5.0, rows*3.8), sharey=sharey)
        axes = axes.ravel() if hasattr(axes, "ravel") else [axes]
        cond_order = cond_order or list(res_by_cond.keys())

        for ax, group in zip(axes, chunked(tags, bonds_per_panel)):
            plotted_any = False
            for bond_tag in group:
                ok = self._overlay_one_bond_on_axis(ax, res_by_cond, bond_tag, to_ns, linewidth, alpha, colors, cond_order, label_with_cond_only)
                plotted_any = plotted_any or ok
            ax.set_xlabel("Lifetime (ns)" if to_ns else "Lifetime (ps)")
            ax.set_ylabel("S(t)")
            ax.set_ylim(0.0, 1.05); ax.grid(True, alpha=0.25); ax.margins(x=0)
            if plotted_any:
                handles, labels = ax.get_legend_handles_labels()
                if len(labels) > 8: handles, labels = handles[:8], labels[:8]
                ax.legend(handles, labels, frameon=False, loc=legend_loc, fontsize=8)
            else:
                ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes); ax.axis("off")

        # hide unused axes
        used = int(np.ceil(len(tags) / bonds_per_panel))
        for ax in axes[used:]:
            ax.axis("off")

        if suptitle: fig.suptitle(suptitle, fontsize=14, y=0.995)
        fig.tight_layout(); plt.show()

    # ---- visualization: occupancy ----
    def plot_occupancy_bars(
        self, res_by_cond: Dict[str, dict], bond_tags: List[str],
        colors: Optional[Dict[str,str]] = None, rotate_xticks: int = 45,
        figsize=(10,6), title: str = "H-bond occupancy (%)"
    ):
        df = self.table_occupancy(res_by_cond, as_percent=True).loc[bond_tags]
        ax = df.plot(kind="bar", figsize=figsize,
                     color=None if colors is None else [colors.get(c, None) for c in df.columns])
        ax.set_ylabel("Occupancy (%)"); ax.set_xlabel(""); ax.set_title(title)
        ax.set_ylim(0, 100); ax.grid(axis="y", alpha=0.2)
        plt.xticks(rotation=rotate_xticks, ha="right"); plt.tight_layout(); plt.show()

    def plot_occupancy_heatmap(
        self, res_by_cond: Dict[str, dict], top_n: int = 50, pair: Optional[Tuple[str,str]] = None,
        title: str = "H-bond occupancy heatmap (%)"
    ):
        df = self.table_occupancy(res_by_cond, as_percent=True)
        if pair is not None:
            a, b = pair; delta = (df.get(a, 0) - df.get(b, 0)).abs()
        else:
            delta = df.max(axis=1) - df.min(axis=1)
        top = delta.sort_values(ascending=False).head(top_n).index
        df_top = df.loc[top]
        plt.figure(figsize=(10, max(4, int(top_n/3))))
        plt.imshow(df_top.values, aspect="auto", interpolation="nearest")
        plt.colorbar(label="Occupancy (%)")
        plt.yticks(range(len(df_top.index)), df_top.index)
        plt.xticks(range(len(df_top.columns)), df_top.columns, rotation=45, ha="right")
        plt.title(title); plt.tight_layout(); plt.show()
