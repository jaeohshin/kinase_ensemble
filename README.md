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
│ │ ├── backbone_*.pdb # Extracted backbone-only structures
│ │ └── gly/ # Mutated glycine-only structures
│ ├── em/ # Energy minimized glycine structures
│ ├── final_backbones/ # Real sequence restored
│ ├── flowpacker/ # Sidechains repacked using FlowPacker
│ └── nvt/ # NVT MD-refined final structures
├── configs/
│ ├── flowpacker_template.yaml # Template config for FlowPacker
│ └── flowpacker_<kinase>.yaml # Auto-generated configs per kinase
├── scripts/
│ ├── run_bioemu.py
│ ├── extract_backbone.py
│ ├── add_glycine_sidechains.py
│ ├── correct_resnames.py
│ ├── minimize_backbones.sh
│ ├── avg_mpscore.sh
│ ├── sampler_pdb.py
│ ├── gromacs_scripts.sh
│ └── run_pipeline.py # Optional: orchestrator script
├── logs/ # Log files for each kinase + step
└── README.md
```


---

## 🧬 Pipeline Overview (with Commands)

Each kinase ensemble is generated through the following steps:

---

### 1. 🧠 BioEmu – Backbone Sampling
Generate backbone-only structures from FASTA.

```bash
python scripts/run_bioemu.py abl1
   → `output/<kinase>/bioemu/samples.xtc`
````
### 2. 🦴 Extract Backbone Atoms

Extract C/N/CA/O atoms from trajectory frames.

````
python scripts/extract_backbone.py output/abl1/bioemu
   → `output/<kinase>/bioemu/backbone_*.pdb`
````
### 3. ⚗️ Mutate All Residues to Glycine

Replace sidechains with glycine and renumber residues.

````
python scripts/add_glycine_sidechains.py output/abl1/bioemu
   → `output/<kinase>/bioemu/gly/backbone_XXX_gly.pdb`
````

### 4. 🔧 Energy Minimization (GROMACS)

Minimize glycine-only structures in vacuum.

````
bash scripts/minimize_backbones.sh output/abl1/bioemu/gly/
   → `output/<kinase>/bioemu/gly/minimized_backbone_XXX_gly.pdb`
````
### 5. 📊 MolProbity Scoring (Pre-repacking)

Evaluate structural quality using MolProbity.

cd /home/deepfold/users/jaeohshin/tools/
bash avg_mpscore.sh
   → Example: `avg 93.96 5.48 0.55 0.46 0.95 265.98`

### 6. 🧬 Restore Real Sequence

Restore the original sequence using residue renaming.

````
python scripts/correct_resnames.py output/abl1/bioemu/gly/
   → `output/<kinase>/final_backbones/backbone_XXX.pdb`
````

### 7. 🎲 FlowPacker – Sidechain Repacking

Run generative sidechain prediction.

    Edit configs/flowpacker_abl1.yaml or generate from a template.
````
python scripts/sampler_pdb.py base abl1

→ Output: tools/flowpacker/samples/abl1/run_1/final_XXX.pdb

→ Sample MolProbity Output: avg 92.16 6.85 0.99 31.78 2.50 265.98
````
### 8. 🌡️ NVT Molecular Dynamics (GROMACS)

Final structure refinement using restrained MD.
````
bash scripts/gromacs_scripts.sh abl1

→ Output: output/abl1/nvt/receptor_XXX.pdb

→ Final MolProbity: avg 92.24 7.15 0.61 1.39 1.35 266.00
````
---

## ✅ Requirements

- Python ≥ 3.10
- [BioEmu](https://github.com/microsoft/bioemu)
- [FlowPacker](https://gitlab.com/mjslee0921/flowpacker)
- [GROMACS](https://www.gromacs.org/) ≥ 2021
- [PDBFixer](https://github.com/openmm/pdbfixer) (OpenMM toolkit)
- MolProbity tools (for quality assessment)
- Conda environments for isolation (recommended)

---

## 🧪 Notes

- Glycine mutation allows for unbiased minimization before sequence-specific packing.
- MolProbity scoring is run **before** and **after** FlowPacker to evaluate clash reduction.
- Structures from the NVT stage are suitable for ligand docking and ensemble-based virtual screening.

