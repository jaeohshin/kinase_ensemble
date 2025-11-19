#!/bin/bash
#SBATCH --job-name=nvt_md_all
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-task=1
#SBATCH --time=96:00:00
#SBATCH --output=logs/slurm_%x_%j.out
#SBATCH --error=logs/slurm_%x_%j.err
#SBATCH --gres-flags=enforce-binding

source ~/miniconda3/etc/profile.d/conda.sh
conda activate gro

KINASES=("igf1r" "jak2" "kpcb" "kit" "lck" "mapk2" "met" "mk01" "mk10")

for kinase in "${KINASES[@]}"; do
    echo "🚀 Starting GROMACS refinement for $kinase at $(date)"
    bash scripts/gromacs_scripts_new.sh "$kinase" 50
    echo "✅ [$kinase] refinement complete at $(date)"
done
