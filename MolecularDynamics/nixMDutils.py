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

