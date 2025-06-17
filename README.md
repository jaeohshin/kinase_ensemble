# Kinase Ensemble Generation Pipeline

This repository contains a modular pipeline to generate ensembles of kinase protein structures, designed for downstream docking and virtual screening experiments.

---

## 📁 Directory Structure
```
kinase_ensemble/
├── input/
│ ├── sequences/ # FASTA files for each kinase
│ └── kinase_list.txt # List of kinase IDs to process
├── output/ # Output directories per kinase
│ └── <kinase>/
│ ├── bioemu/ # Raw backbone samples from BioEmu
│ ├── pdbfixer/ # Sidechains added via PDBFixer
│ ├── em/ # Energy minimized structures (GROMACS)
│ ├── flowpacker/ # Sidechains repacked using FlowPacker
│ └── nvt/ # NVT MD-refined final structures
├── configs/
│ ├── flowpacker_template.yaml # Template config for FlowPacker
│ └── flowpacker_<kinase>.yaml # Auto-generated configs per kinase
├── scripts/
│ ├── run_bioemu.py
│ ├── run_pdbfixer.py
│ ├── run_em.py
│ ├── run_flowpacker.py
│ ├── run_nvt.sh
│ ├── run_pipeline.py # Orchestrator for full pipeline
│ └── generate_flowpacker_yaml.py # Optional: templating utility
├── logs/ # Log files for each kinase + step
└── README.md
```

---

## 🧬 Pipeline Overview

Each kinase structure is generated from a FASTA sequence through five key stages:

1. **BioEmu**: Sample backbone-only structures.
2. **PDBFixer**: Add initial sidechains for minimization.
3. **Energy Minimization**: Refine structures using GROMACS.
4. **FlowPacker**: Repack sidechains using a generative model.
5. **NVT MD**: Final refinement via restrained molecular dynamics.

## ✅ Requirements

- Python ≥ 3.10
- BioEmu, FlowPacker installed via pip or local source
- GROMACS ≥ 2021
- PDBFixer (from OpenMM toolkit)
- conda environments per step recommended
