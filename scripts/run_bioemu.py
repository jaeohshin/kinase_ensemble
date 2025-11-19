import os
import sys
import subprocess

template = """#!/bin/bash
# Direct run script for BioEmu (no SLURM)

# Load conda and activate env
source ~/miniconda3/etc/profile.d/conda.sh
conda activate bioemu

# Extract sequence from FASTA
sequence=$(grep -v "^>" {fasta} | tr -d '\\n')

# Run BioEmu
python -m bioemu.sample \\
    --sequence "$sequence" \\
    --num_samples 50 \\
    --output_dir {out_dir} \\
    --cache_embeds_dir /store/jaeohshin/.bioemu_embeds_cache \\
    --filter_samples False
"""

def make_local_script(kinase: str):
    fasta = f"input/fasta/{kinase}.fasta"
    out_dir = f"output/{kinase}/bioemu"
    log_dir = f"logs/{kinase}"

    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)

    if not os.path.isfile(fasta):
        print(f"[ERROR] FASTA not found: {fasta}")
        return

    script_path = os.path.join(out_dir, "run_bioemu.sh")
    with open(script_path, "w") as f:
        f.write(template.format(kinase=kinase, fasta=fasta, out_dir=out_dir))

    os.chmod(script_path, 0o755)
    print(f"[INFO] Execution script created at {script_path}")
    return script_path

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python scripts/run_bioemu.py <kinase>")
        sys.exit(1)

    kinase = sys.argv[1].lower()
    script_path = make_local_script(kinase)

    if script_path:
        print(f"[INFO] Running BioEmu for {kinase}")
        subprocess.run([script_path])
