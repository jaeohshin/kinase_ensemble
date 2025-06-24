import sys
import yaml
from pathlib import Path

TEMPLATE = {
    'mode': 'vf',
    'data': {
        'data': 'bc40',
        'train_path': '',
        'cluster_path': '',
        'test_path': '',  # will be filled dynamically
        'min_length': 40,
        'max_length': 512,
        'edge_type': 'knn',
        'max_radius': 16.0,
        'max_neighbors': 30
    },
    'ckpt': '/store/jaeohshin/tools/flowpacker/checkpoints/cluster.pth',
    'conf_ckpt': None,  # can be replaced with real path if used
    'sample': {
        'batch_size': 1,
        'n_samples': 1,
        'use_ema': True,
        'eps': 2.0e-3,
        'save_trajectory': False,
        'coeff': 5.0,
        'num_steps': 10
    },
    'output_dir': ''  # will be filled dynamically
}

CONFIG_DIR = "../configs"

def generate_yaml(kinase):
    config = TEMPLATE.copy()
    config['data'] = TEMPLATE['data'].copy()
    config['sample'] = TEMPLATE['sample'].copy()

    config['data']['test_path'] = f"/store/jaeohshin/work/kine/output/{kinase}/em"
    config['output_dir'] = f"/store/jaeohshin/work/kine/output/{kinase}/flowpacker"

    config_path = Path(CONFIG_DIR) / f"flowpacker_{kinase}.yaml"
    config_path.parent.mkdir(parents=True, exist_ok=True)

    with open(config_path, 'w') as f:
        yaml.dump(config, f, sort_keys=False, default_flow_style=False)

    print(f"[✓] YAML generated: {config_path}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python generate_flowpacker_yaml.py <kinase>")
        sys.exit(1)

    kinase = sys.argv[1]
    generate_yaml(kinase)

