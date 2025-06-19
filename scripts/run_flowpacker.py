import argparse
import os
import shutil
import subprocess
from pathlib import Path

def run_flowpacker(kinase):
    project_root = Path(__file__).resolve().parents[1]

    input_dir = project_root / "output" / kinase / "bioemu" / "final_backbones"
    flowpacker_dir = project_root / "output" / kinase / "flowpacker"
    config_template = project_root / "configs" / "flowpacker_template.yaml"
    config_path = project_root / "configs" / f"flowpacker_{kinase}.yaml"
    sample_output_dir = Path("/store/jaeohshin/tools/flowpacker/samples") / f"{kinase}_flowpacker" / "run_1"

    flowpacker_dir.mkdir(parents=True, exist_ok=True)
    flowpacker_input_dir = sample_output_dir.parent
    flowpacker_input_dir.mkdir(parents=True, exist_ok=True)

    # Copy corrected PDBs to FlowPacker input directory
    print(f"[INFO] Copying PDBs from {input_dir} to {flowpacker_input_dir}...")
    for pdb_file in sorted(input_dir.glob("minimized_backbone_*_corrected.pdb")):
        shutil.copy(pdb_file, flowpacker_input_dir / pdb_file.name)

    # Generate config from template
    with open(config_template) as fin:
        config_text = fin.read()
        config_text = config_text.replace("REPLACE_ME_WITH_TEST_PATH", str(flowpacker_input_dir))
        config_text = config_text.replace("./checkpoints/cluster.pth", "/store/jaeohshin/tools/flowpacker/checkpoints/cluster.pth")

    # Write config to FlowPacker expected location
    flowpacker_config_dir = Path("/store/jaeohshin/tools/flowpacker/config/inference")
    flowpacker_config_dir.mkdir(parents=True, exist_ok=True)
    final_config_path = flowpacker_config_dir / "base.yaml"
    with open(final_config_path, "w") as fout:
        fout.write(config_text)

    # Run FlowPacker
    cmd = [
        "python", "/store/jaeohshin/tools/flowpacker/sampler_pdb.py",
        "base", kinase,
        "--use_gt_masks", "True"
    ]

    print("[INFO] Running FlowPacker...")
    subprocess.run(cmd, check=True, cwd="/store/jaeohshin/tools/flowpacker")

    # Move outputs to flowpacker dir
    for pdb_file in sample_output_dir.glob("*.pdb"):
        shutil.move(str(pdb_file), flowpacker_dir / pdb_file.name)

    print(f"[INFO] FlowPacker run complete. Outputs saved to: {flowpacker_dir}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("kinase", help="Kinase name (e.g., abl1)")
    args = parser.parse_args()
    run_flowpacker(args.kinase)
