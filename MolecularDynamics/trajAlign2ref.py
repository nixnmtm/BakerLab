'''
@author:     Nixon Raj

@copyright:  2023 Nationwide Childrens Hospital, Abigail Wexner Research Institute

@license:    GPL v. 3

@contact:    nixon.raj@nationwidechildrens.org
@deffield    updated: Updated
'''

import argparse
import os
import sys

import MDAnalysis as mda
import MDAnalysis.analysis.align as mdalign
import os.path as path

program_version_message="v2.0"

def parse_args():
    """Parse command-line arguments
    """
    parser = argparse.ArgumentParser(description="Align Trajectories",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-V", "--version", action="version",
                        version=program_version_message)
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Write additional information to screen")

    # I/O group
    iogrp = parser.add_argument_group(title="File I/O")
    iogrp.add_argument("-s", dest="structin",
                       help="Structure file, input [default: %(default)s]",
                       default="conf.gro", metavar="PSF")
    iogrp.add_argument("-f", dest="trajin",
                       help="Trajectory file, input [default: %(default)s]",
                       default="traj.xtc", metavar="TRJ")
    iogrp.add_argument("-o", dest="trajout",
                       help="Trajectory file, output [default: %(default)s]",
                       default=path.join(os.curdir, "traj.aligned.xtc"),
                       metavar="TRJ")
    iogrp.add_argument("-rs", "--refS", dest="ref_struc", metavar="REF_STRUC",
                       help="Reference structure for alignment PSF")
    iogrp.add_argument("-rc", "--refC", dest="ref_coor", metavar="REF_COOR",
                       help="Reference structure for alignment CRD")
    iogrp.add_argument("-sel", "--select", dest="selection", type=str, metavar="SELECTION",
                       help="example name Ca and segid 3CM5")
    # Process arguments
    args = parser.parse_args()

    print(parser.description.replace("\nUSAGE", ""))
    print(args, "\n")

    return args

def align_traj(args):
    """Merge coordinate files into a single trajectory file.

    :param args: command-line arguments
    """
    # Setup reference structure for alignment
    if args.ref_struc and args.ref_coor: # set the given reference
        print("Using the given reference structure")
        refpsf = args.ref_struc
        reftraj = args.ref_coor
    else:
        refpsf = args.structin
        reftraj = args.trajin

    ref_univ = mda.Universe(refpsf, reftraj)

    # Create universe consisting of all coordinates.
    test_univ = mda.Universe(args.structin, args.trajin)

    # Align the coordinates to the reference structure and write to a
    # trajectory file.

    select = "all"
    if args.selection:
        select = args.selection
    if args.verbose:
        aligned = mdalign.AlignTraj(test_univ, ref_univ, select=select,
                filename=args.trajout, verbose=True)
    else:
        aligned = mdalign.AlignTraj(test_univ, ref_univ, select=select,
                filename=args.trajout)
    aligned.run()

def main():  # IGNORE:C0111
    """Command line options."""
    try:
        # Setup argument parser
        args = parse_args()
        align_traj(args)
        return os.EX_OK
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return "Error"

if __name__ == "__main__":
    main()