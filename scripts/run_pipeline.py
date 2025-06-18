# scripts/run_pipeline.py
import os
import subprocess

KINASE_LIST = "input/kinase_list.txt"

def run_cmd(cmd):
    print(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)

with open(KINASE_LIST) as f:
    kinases = [line.strip() for line in f if line.strip()]

for kinase in kinases:
    print(f"\n=== Processing {kinase} ===")
    
    # Step 1: BioEmu
    run_cmd(f"python scripts/run_bioemu.py {kinase}")
    
    # Step 2: EM
    run_cmd(f"python scripts/run_em.py {kinase}")
    
    # Step 3: FlowPacker
    run_cmd(f"python scripts/run_flowpacker.py {kinase}")
    
    # Step 4: NVT
    run_cmd(f"bash scripts/run_nvt.md.sh {kinase}")
