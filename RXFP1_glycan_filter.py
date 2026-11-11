"""Glycan-proximity filter for de novo RXFP1 binder designs.

Excludes designs whose binder chain comes within a cutoff distance of the
glycan N-linked to N335 (N303 in signal-peptide-subtracted RXFP1 numbering).

The glycan is represented by a single pseudo-atom, stored as a water (HOH)
residue in a marker PDB. Its position was assigned from density associated
with N335 in the low-resolution cryo-EM map of RXFP1 (EMD-26004).

The only required argument is the cutoff in angstroms. Designs whose
closest binder atom sits at or beyond the cutoff are copied to
filtered_pdbs_<cutoff>A/; the rest go to failed_pdbs_<cutoff>A/. Inputs are
never modified. A per-design distance table and the pass count are printed
to stdout, so redirect it if you want a record of the triage.

Usage:
    # Defaults: 6.5 A cutoff, binder in chain B, all atoms, reading
    # ./*/Accepted/*.pdb (the BindCraft output layout) and writing the two
    # output folders into the current directory.
    python rxfp1_glycan_filter.py 6.5

    # Same, but for designs whose binder is chain A rather than chain B.
    python rxfp1_glycan_filter.py 6.5 --binder-chain A

    # Point at a different set of designs, and measure backbone atoms only
    # (N, CA, C, O), ignoring sidechains. Quote the glob so the shell passes
    # it through rather than expanding it.
    python rxfp1_glycan_filter.py 6.5 --input "./designs/*.pdb" --backbone-only

    # Marker file elsewhere, outputs written somewhere other than the cwd.
    python rxfp1_glycan_filter.py 6.5 --marker ../RXFP1_Glyco_marker.pdb \
        --outdir ./triage

Run with --help for the full argument list and defaults.

Note the binder chain differs between design pipelines: BindCraft and
RFdiffusion do not use the same convention, so check your inputs before
relying on the default. A design with no atoms in the named chain is
skipped with a warning rather than silently passing the filter.
"""

import argparse
import glob
import os
import shutil
import sys

import numpy as np
from Bio import PDB

BACKBONE_ATOMS = {"N", "CA", "C", "O"}


def get_marker_coordinate(marker_pdb):
    """Return the coordinate of the pseudo-atom marking the glycan position.

    The marker is the first atom of the first HOH residue in the file.
    """
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("marker", marker_pdb)

    waters = [
        residue
        for residue in structure.get_residues()
        if residue.get_resname() == "HOH"
    ]

    if not waters:
        raise ValueError(f"No HOH residue found in marker file {marker_pdb}")
    if len(waters) > 1:
        print(
            f"Warning: {len(waters)} HOH residues in {marker_pdb}; "
            "using the first.",
            file=sys.stderr,
        )

    return next(waters[0].get_atoms()).coord


def min_distance_to_marker(pdb_file, marker_coord, binder_chain, backbone_only=False):
    """Minimum distance from any binder-chain atom to the glycan pseudo-atom."""
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("design", pdb_file)

    target = binder_chain.strip().upper()
    min_dist = float("inf")

    for chain in structure.get_chains():
        if chain.id.strip().upper() != target:
            continue
        for atom in chain.get_atoms():
            if backbone_only and atom.get_name() not in BACKBONE_ATOMS:
                continue
            dist = float(np.linalg.norm(atom.coord - marker_coord))
            min_dist = min(min_dist, dist)

    if min_dist == float("inf"):
        print(
            f"Warning: no chain '{binder_chain}' atoms in {pdb_file}; skipping.",
            file=sys.stderr,
        )

    return min_dist


def filter_designs(
    input_pattern,
    marker_pdb,
    passed_dir,
    failed_dir,
    threshold,
    binder_chain,
    backbone_only=False,
):
    marker_coord = get_marker_coordinate(marker_pdb)
    print(f"Glycan pseudo-atom at {marker_coord}")

    pdb_files = sorted(glob.glob(input_pattern))
    if not pdb_files:
        print(f"No PDB files matched {input_pattern}", file=sys.stderr)
        return []

    print(f"Found {len(pdb_files)} design(s) matching {input_pattern}")

    os.makedirs(passed_dir, exist_ok=True)
    os.makedirs(failed_dir, exist_ok=True)

    passed = []

    for pdb_file in pdb_files:
        min_dist = min_distance_to_marker(
            pdb_file, marker_coord, binder_chain, backbone_only
        )
        if min_dist == float("inf"):
            continue

        keep = min_dist >= threshold
        dest_dir = passed_dir if keep else failed_dir
        dest_file = os.path.join(dest_dir, os.path.basename(pdb_file))

        if os.path.abspath(pdb_file) != os.path.abspath(dest_file):
            shutil.copy(pdb_file, dest_file)

        print(f"{os.path.basename(pdb_file)}\t{min_dist:.2f} A\t"
              f"{'PASS' if keep else 'FAIL'}")

        if keep:
            passed.append(dest_file)

    print(f"\n{len(passed)}/{len(pdb_files)} passed at {threshold} A")
    print(f"Passed -> {passed_dir}")
    print(f"Failed -> {failed_dir}")

    return passed


def main():
    parser = argparse.ArgumentParser(
        description="Exclude binder designs that clash with the RXFP1 N335 glycan."
    )
    parser.add_argument(
        "threshold", type=float, help="Minimum allowed distance in angstroms (e.g. 6.5)"
    )
    parser.add_argument(
        "--input",
        default="./*/Accepted/*.pdb",
        help="Glob pattern for design PDBs (default: ./*/Accepted/*.pdb)",
    )
    parser.add_argument(
        "--marker",
        default="RXFP1_Glyco_marker.pdb",
        help="Marker PDB containing the glycan pseudo-atom as a HOH residue",
    )
    parser.add_argument(
        "--binder-chain",
        default="B",
        help="Chain ID of the binder (default: B; differs between pipelines)",
    )
    parser.add_argument(
        "--backbone-only",
        action="store_true",
        help="Measure backbone atoms only (default: all atoms)",
    )
    parser.add_argument(
        "--outdir",
        default=".",
        help="Directory in which to create the passed/failed folders",
    )
    args = parser.parse_args()

    tag = f"{args.threshold:.1f}A"
    passed_dir = os.path.join(args.outdir, f"filtered_pdbs_{tag}")
    failed_dir = os.path.join(args.outdir, f"failed_pdbs_{tag}")

    filter_designs(
        input_pattern=args.input,
        marker_pdb=args.marker,
        passed_dir=passed_dir,
        failed_dir=failed_dir,
        threshold=args.threshold,
        binder_chain=args.binder_chain,
        backbone_only=args.backbone_only,
    )


if __name__ == "__main__":
    main()
