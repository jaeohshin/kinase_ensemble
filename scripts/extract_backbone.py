#!/usr/bin/env python3

"""
Extracts backbone structures from a BioEmu-generated trajectory
and saves them as individual PDB files in the same directory.

Usage:
    python scripts/extract_backbones.py output/abl1/bioemu
"""

import os
import sys
import mdtraj as md

def extract_backbones(bioemu_dir):
    top = os.path.join(bioemu_dir, "topology.pdb")
    xtc = os.path.join(bioemu_dir, "samples.xtc")

    if not os.path.isfile(top) or not os.path.isfile(xtc):
        print(f"[ERROR] Missing topology or trajectory in: {bioemu_dir}")
        return

    print(f"[INFO] Loading trajectory from {xtc}")
    traj = md.load(xtc, top=top)

    for i, frame in enumerate(traj):
        out_path = os.path.join(bioemu_dir, f"backbone_{i:03d}.pdb")
        frame.save_pdb(out_path)
        print(f"[INFO] Saved: {out_path}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python scripts/extract_backbones.py <bioemu_output_dir>")
        sys.exit(1)

    extract_backbones(sys.argv[1])

