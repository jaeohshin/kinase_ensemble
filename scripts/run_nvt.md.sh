#!/bin/bash
# scripts/run_nvt.md.sh

kinase=$1
flow_dir="output/${kinase}/flowpacker"
cd $flow_dir

# Assume final PDB is named like flowpacked_${kinase}.pdb
pdb="flowpacked_${kinase}.pdb"
gmx pdb2gmx -f $pdb -o processed.gro -water tip3p
gmx editconf -f processed.gro -o newbox.gro -c -d 1.0 -bt cubic
gmx solvate -cp newbox.gro -cs spc216.gro -o solv.gro -p topol.top
gmx grompp -f nvt.mdp -c solv.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt
