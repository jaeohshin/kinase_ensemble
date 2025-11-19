#!/usr/bin/env python3
import os
import sys
import subprocess
import glob
import shutil

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
    --num_samples {num_samples} \\
    --output_dir {temp_dir} \\
    --cache_embeds_dir /store/jaeohshin/.bioemu_embeds_cache \\
    --filter_samples False
"""

def get_missing_numbers(out_dir: str, target_total: int = 50) -> tuple:
    """Find which sample numbers (000-049) are missing."""
    existing_nums = set()
    
    # Check for all possible PDB files
    for pdb_file in glob.glob(os.path.join(out_dir, "*.pdb")):
        basename = os.path.basename(pdb_file)
        # Try to extract 3-digit number from filename
        for i in range(target_total):
            num_str = f"{i:03d}"
            if num_str in basename:
                existing_nums.add(i)
                break
    
    missing = [i for i in range(target_total) if i not in existing_nums]
    existing = sorted(existing_nums)
    
    return missing, existing

def merge_new_structures(temp_dir: str, out_dir: str, missing_nums: list):
    """Copy and rename new structures from temp to fill missing numbers."""
    new_files = sorted(glob.glob(os.path.join(temp_dir, "*.pdb")))
    
    if not new_files:
        print(f"[ERROR] No PDB files generated in {temp_dir}")
        return False
    
    print(f"[INFO] Found {len(new_files)} newly generated structures")
    
    # Map new files to missing numbers
    for new_file, target_num in zip(new_files, missing_nums):
        # Get the original filename pattern to preserve it
        basename = os.path.basename(new_file)
        
        # Try to detect the naming pattern and create target filename
        # Common patterns: sample_000.pdb, 000.pdb, relaxed_000.pdb, etc.
        if "sample_" in basename:
            target_file = os.path.join(out_dir, f"sample_{target_num:03d}.pdb")
        elif basename.startswith("relaxed_"):
            target_file = os.path.join(out_dir, f"relaxed_{target_num:03d}.pdb")
        else:
            # Default to sample_ prefix
            target_file = os.path.join(out_dir, f"sample_{target_num:03d}.pdb")
        
        shutil.copy2(new_file, target_file)
        print(f"[INFO] Copied {os.path.basename(new_file)} -> {os.path.basename(target_file)}")
    
    # Clean up temp directory
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
        print(f"[INFO] Cleaned up temporary directory")
    
    return True

def make_local_script(kinase: str, target_total: int = 50):
    fasta = f"input/fasta/{kinase}.fasta"
    out_dir = f"output/{kinase}/bioemu"
    temp_dir = f"output/{kinase}/bioemu_temp"
    log_dir = f"logs/{kinase}"
    
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)
    
    if not os.path.isfile(fasta):
        print(f"[ERROR] FASTA not found: {fasta}")
        return None
    
    # Check what's missing
    missing, existing = get_missing_numbers(out_dir, target_total)
    
    print(f"[INFO] Existing structures: {len(existing)}/{target_total}")
    
    if not missing:
        print(f"[INFO] All {target_total} structures already exist. Skipping.")
        return None
    
    print(f"[INFO] Missing {len(missing)} structures: {missing[:5]}" + 
          (f"...{missing[-1]}" if len(missing) > 5 else ""))
    
    # Create temp directory for new generation
    os.makedirs(temp_dir, exist_ok=True)
    
    script_path = os.path.join(temp_dir, "run_bioemu.sh")
    with open(script_path, "w") as f:
        f.write(template.format(
            kinase=kinase,
            fasta=fasta,
            temp_dir=temp_dir,
            num_samples=len(missing)
        ))
    
    os.chmod(script_path, 0o755)
    print(f"[INFO] Execution script created at {script_path}")
    
    return script_path, temp_dir, out_dir, missing

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python scripts/run_bioemu.py <kinase>")
        sys.exit(1)
    
    kinase = sys.argv[1].lower()
    result = make_local_script(kinase)
    
    if result is None:
        print(f"[INFO] No action needed for {kinase}")
        sys.exit(0)
    
    script_path, temp_dir, out_dir, missing = result
    
    print(f"[INFO] Running BioEmu for {kinase} (generating {len(missing)} structures)")
    proc_result = subprocess.run([script_path])
    
    if proc_result.returncode == 0:
        print(f"[INFO] BioEmu completed successfully")
        success = merge_new_structures(temp_dir, out_dir, missing)
        if success:
            print(f"[SUCCESS] ✓ {kinase} now has 50 structures")
        else:
            print(f"[ERROR] Failed to merge structures for {kinase}")
            sys.exit(1)
    else:
        print(f"[ERROR] BioEmu failed for {kinase}")
        sys.exit(1)
