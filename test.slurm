import os
import sys

# Template for the SLURM job script
template = """#!/bin/bash
#SBATCH --job-name=bioemu_{kinase}
#SBATCH --output={log_dir}/bioemu_%j.out
#SBATCH --error={log_dir}/bioemu_%j.err
#SBATCH --time=02:00:00
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

# Load conda and activate env
source ~/miniconda3/etc/profile.d/conda.sh
conda activate bioemu

# Set temp dir
export TMPDIR=/store/jaeohshin/tmp_bioemu/$SLURM_JOB_ID
mkdir -p $TMPDIR

# Extract sequence from FASTA
sequence=$(grep -v "^>" {fasta} | tr -d '\\n')

# Run BioEmu
python -m bioemu.sample \\
    --sequence "$sequence" \\
    --num_samples 2 \\
    --output_dir {out_dir} \\
    --cache_embeds_dir /store/jaeohshin/.bioemu_embeds_cache \\
    --filter_samples False
"""

def make_slurm_script(kinase: str):
    fasta = f"input/fasta/{kinase}_trimmed.fasta"
    out_dir = f"output/{kinase}/bioemu"
    log_dir = f"logs/{kinase}"

    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)

    if not os.path.isfile(fasta):
        print(f"[ERROR] FASTA not found: {fasta}")
        return

    job_script_path = os.path.join(out_dir, "bioemu_job.sh")
    with open(job_script_path, "w") as f:
        f.write(template.format(kinase=kinase, fasta=fasta, out_dir=out_dir, log_dir=log_dir))

    print(f"[INFO] SLURM script created at {job_script_path}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python scripts/run_bioemu.py <kinase>")
        sys.exit(1)

    kinase = sys.argv[1].lower()
    make_slurm_script(kinase)

