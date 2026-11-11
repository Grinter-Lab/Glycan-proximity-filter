# rxfp1-glycan-filter

Geometric filter excluding de novo RXFP1 binder designs that clash with the N335 glycan.

RXFP1 carries an N-linked glycan at N335 (N303 in the signal-peptide-subtracted numbering scheme) close to the H2 relaxin binding site. Design pipelines model the receptor without glycans, so designs can score well while occupying space the glycan takes up. This script removes them.

## Method

The glycan is represented by a single pseudo-atom rather than a modelled oligosaccharide. Its position was assigned from density associated with N335 in the low-resolution cryo-EM map of RXFP1 ([EMD-26004](https://www.ebi.ac.uk/emdb/EMD-26004)), and is stored as a `HOH` residue in `RXFP1_Glyco_marker.pdb`.

For each design, the script computes the distance from every atom of the binder chain to that pseudo-atom and keeps the minimum. Designs are retained if that minimum is at or above the cutoff. A 6.5 Å cutoff approximates the excluded volume of the proximal N-acetylglucosamine plus a clash buffer.

## Usage

```
python rxfp1_glycan_filter.py 6.5
```

Passing designs are copied to `filtered_pdbs_6.5A/`, the rest to `failed_pdbs_6.5A/`. Inputs are never modified. A per-design distance table and the pass count are printed to stdout.

| Argument          | Default                  | Description                        |
| ----------------- | ------------------------ | ---------------------------------- |
| `threshold`       | —                        | Cutoff in angstroms (required)     |
| `--input`         | `./*/Accepted/*.pdb`     | Glob for design PDBs               |
| `--marker`        | `RXFP1_Glyco_marker.pdb` | Marker PDB holding the pseudo-atom |
| `--binder-chain`  | `B`                      | Chain ID of the binder             |
| `--backbone-only` | off                      | Measure N, CA, C, O only           |
| `--outdir`        | `.`                      | Where to create the output folders |

**Check the binder chain.** BindCraft and RFdiffusion outputs do not use the same chain convention, so confirm which chain holds the binder before running. A design with no atoms in the named chain is skipped with a warning rather than silently passing.

## Requirements

Python 3, [Biopython](https://biopython.org/), NumPy.

```
pip install biopython numpy
```

## Citation

If you use this script, please cite the Zenodo archive (DOI 10.5281/zenodo.24126030) and the paper it was written for. 

A complete RXFP1–relaxin interaction model unlocks the design of potent mini-protein modulators
Janik Clement, Tim Lkhagvajargal, Bradley L. Hoare, Tiffany Myint, Daniel R. Fox, Chunxiao Wang, Gavin J. Knott, Ross A.D. Bathgate, Rhys Grinter
bioRxiv 2026.06.19.733483; doi: https://doi.org/10.64898/2026.06.19.733483
