from itertools import groupby

def generate_gromacs_ndx_selection_from_residues(residue_set):
    """
    Takes a set of integer residue numbers and returns a GROMACS-style
    make_ndx selection string with ranges grouped for compactness.
    """
    residues = sorted(residue_set)
    ranges = []
    for k, g in groupby(enumerate(residues), lambda x: x[1] - x[0]):
        group = list(g)
        start = group[0][1]
        end = group[-1][1]
        if start == end:
            ranges.append(f"r {start}")
        else:
            ranges.append(f"r {start}-{end}")
    return " | ".join(ranges)

def generate_vmd_selection_from_residues(residues) -> str:
    """
    Build a compact VMD selection string from residue numbers.

    VMD accepts selections like:
      - "resid 5 9 12"                  (discrete residues)
      - "resid 5 to 9 12 20 to 25"     (ranges + discrete)
    Multiple residues/ranges after a single 'resid' are implicitly OR'ed.

    Args:
        residues: Iterable of residue indices (ints). Duplicates are ignored.

    Returns:
        A string like "resid 5 7 to 10 42". If no residues are given,
        returns "none" (a valid empty selection in VMD).

    Raises:
        TypeError: if any residue cannot be interpreted as an int.
    """
    # Normalize & validate
    try:
        uniq: Set[int] = {int(r) for r in residues}
    except (ValueError, TypeError) as e:
        raise TypeError("All residues must be integers or castable to int.") from e

    if not uniq:
        return "none"

    sorted_res = sorted(uniq)

    # Collapse consecutive runs into ranges
    ranges = []
    start = prev = sorted_res[0]
    for r in sorted_res[1:]:
        if r == prev + 1:
            prev = r
            continue
        # Close current run
        ranges.append((start, prev))
        start = prev = r
    ranges.append((start, prev))

    # Format: single numbers as N, runs as "A to B"
    parts = []
    for a, b in ranges:
        if a == b:
            parts.append(str(a))
        else:
            parts.append(f"{a} to {b}")

    return "resid " + " ".join(parts)


    
