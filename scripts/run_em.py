# scripts/run_em.py
import sys
import os
import subprocess

kinase = sys.argv[1]
input_dir = f"output/{kinase}/bioemu"
output_dir = f"output/{kinase}/em"
os.makedirs(output_dir, exist_ok=True)

# Assume you apply PDBFixer first
for pdb in os.listdir(input_dir):
    if pdb.endswith(".pdb"):
        in_path = os.path.join(input_dir, pdb)
        fixed_path = os.path.join(output_dir, f"fixed_{pdb}")
        
        fixer_cmd = f"python scripts/fix_pdb.py {in_path} {fixed_path}"  # Make this script if needed
        subprocess.run(fixer_cmd, shell=True, check=True)

        # EM with GROMACS (example — adapt to your .mdp + path setup)
        gro_cmd = f"gmx grompp -f em.mdp -c {fixed_path} -o em.tpr && gmx mdrun -deffnm em"
        subprocess.run(gro_cmd, shell=True, check=True)
