import os
import sys
import shutil
import subprocess
from pathlib import Path

KINASE_CONFIG_DIR = "/store/jaeohshin/work/kine/configs"
FLOWPACKER_ROOT = "/store/jaeohshin/tools/flowpacker"
FLOWPACKER_CONFIG_DIR = f"{FLOWPACKER_ROOT}/config/inference"
FLOWPACKER_SCRIPT = f"{FLOWPACKER_ROOT}/sampler_pdb.py"

def run_flowpacker(kinase):
    config_path = Path(KINASE_CONFIG_DIR) / f"flowpacker_{kinase}.yaml"
    target_path = Path(FLOWPACKER_CONFIG_DIR) / f"{kinase}.yaml"

    if not config_path.exists():
        print(f"[ERROR] Config not found: {config_path}")
        sys.exit(1)

    if not Path(FLOWPACKER_SCRIPT).exists():
        print(f"[ERROR] FlowPacker entry script not found: {FLOWPACKER_SCRIPT}")
        sys.exit(1)

    print(f"[INFO] Copying config to FlowPacker: {target_path}")
    target_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(config_path, target_path)

    print(f"[INFO] Running FlowPacker on kinase: {kinase}")
    env = os.environ.copy()
    env["PYTHONPATH"] = FLOWPACKER_ROOT  # ensure imports work
    
    subprocess.run(
        ["python", "sampler_pdb.py", kinase, "run1"],
        cwd=FLOWPACKER_ROOT,
        env=env,
        check=True
    )


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python run_flowpacker.py <kinase>")
        sys.exit(1)

    run_flowpacker(sys.argv[1])
