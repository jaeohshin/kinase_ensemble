# scripts/run_flowpacker.py
import sys
import os
import subprocess
import shutil

kinase = sys.argv[1]
em_dir = f"output/{kinase}/em"
out_dir = f"output/{kinase}/flowpacker"
os.makedirs(out_dir, exist_ok=True)

pdb_to_pack = os.path.join(em_dir, "em.gro")  # Or whatever the last EM result is

# Prepare config YAML
flow_yaml = f"configs/{kinase}.yaml"
shutil.copy("configs/template.yaml", flow_yaml)

# Modify test_path line in YAML to point to `pdb_to_pack` dir
with open(flow_yaml, 'r') as f:
    lines = f.readlines()
with open(flow_yaml, 'w') as f:
    for line in lines:
        if line.strip().startswith("test_path:"):
            f.write(f"  test_path: {em_dir}\n")
        else:
            f.write(line)

cmd = f"python /path/to/flowpacker/sampler_pdb.py base {kinase} --use_gt_masks True"
subprocess.run(cmd, shell=True, check=True)
