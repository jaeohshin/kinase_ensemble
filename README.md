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

## 🧬 Pipeline Overview

Each kinase ensemble is generated through the following stages:

1. **BioEmu**  
   Sample diverse backbone-only structures from a single FASTA sequence.  
   → `output/<kinase>/bioemu/samples.xtc`

2. **Extract Backbone**  
   Extract C/N/CA/O atoms from each structure frame.  
   → `output/<kinase>/bioemu/backbone_*.pdb`

3. **Mutate to Glycine**  
   All residues are mutated to glycine, and residue indices are renumbered from 1.  
   → `output/<kinase>/bioemu/gly/backbone_XXX_gly.pdb`

4. **Energy Minimization**  
   Structures are minimized using GROMACS in vacuum (non-periodic) settings.  
   → `output/<kinase>/bioemu/gly/minimized_backbone_XXX_gly.pdb`

5. **MolProbity Score Evaluation**  
   Quality of minimized structures is assessed (clash score, rotamer outliers, etc.).  
   → Example: `avg 93.96 5.48 0.55 0.46 0.95 265.98`

6. **Restore Real Sequence**  
   Replace glycine-only residues with the correct sequence from input FASTA.  
   → `output/<kinase>/final_backbones/backbone_XXX.pdb`

7. **FlowPacker**  
   Use a generative model to repack full-atom sidechains.  
   → `tools/flowpacker/samples/<kinase>/run_1/final_XXX.pdb`  
   → MolProbity scores are evaluated again.

8. **NVT Molecular Dynamics**  
   Perform restrained NVT MD simulations for final refinement.  
   → `output/<kinase>/nvt/receptor_XXX.pdb`

---

## ✅ Requirements

- Python ≥ 3.10
- [BioEmu](https://github.com/microsoft/bioemu)
- [FlowPacker](https://github.com/microsoft/flowpacker)
- [GROMACS](https://www.gromacs.org/) ≥ 2021
- [PDBFixer](https://github.com/openmm/pdbfixer) (OpenMM toolkit)
- MolProbity tools (for quality assessment)
- Conda environments for isolation (recommended)

---

## 🧪 Notes

- Glycine mutation allows for unbiased minimization before sequence-specific packing.
- MolProbity scoring is run **before** and **after** FlowPacker to evaluate clash reduction.
- Structures from the NVT stage are suitable for ligand docking and ensemble-based virtual screening.

