#!/bin/bash
set -e
source ~/miniconda3/etc/profile.d/conda.sh

#KINASES=("abl1" "akt1" "akt2" "braf" "cdk2" "csf1r" "egfr" "fak1" "fgfr1" "igf1r" "jak2" "kpcb" "kit" "lck" "mapk2" "met" "mk01" "mk10" "mk14" "mp2k1" "plk1" "rock1" "src" "tgfr1" "vgfr2" "wee1")
KINASES=("abl1" "akt1" "akt2" "braf" "cdk2" "csf1r" "egfr" "fgfr1" "jak2" "kpcb" "kit" "lck" "mapk2" "met" "mk01" "mk10" "mp2k1" "plk1" "rock1" "src" "tgfr1" "vgfr2")


for kinase in "${KINASES[@]}"; do
    echo "🔵 Processing: $kinase"


    # Step 1 – BioEmu
    echo "→ Step 1: BioEmu"
#    conda activate bioemu
#    python scripts/run_bioemu.py "$kinase"

    # Step 2 – Extract backbone
#    echo "→ Step 2: Extract backbone"
#    python scripts/extract_backbone.py "output/$kinase/bioemu"

    # Step 3 – Glycine mutation
#    echo "→ Step 3: Add Gly"
#    python scripts/add_glycine_sidechains.py "output/$kinase/bioemu"

    # Step 4 – Energy Minimization
    echo "→ Step 4: Minimization"
#    conda activate gro
#    bash scripts/minimize_backbones.sh "output/$kinase/bioemu/gly/"

    # Step 5 – MolProbity scoring (optional, insert here if needed)

    # Step 6 – Restore real sequence
    echo "→ Step 6: Restore sequence"
#    python scripts/correct_resnames.py "$kinase"


    # Step 7 – FlowPacker
    echo "→ Step 7: FlowPacker"
    conda activate flowpacker
    FLOWPACKER_DIR="/store/jaeohshin/tools/flowpacker"
    BASE_CONFIG="${FLOWPACKER_DIR}/config/inference/base.yaml"
    KINASE_CONFIG="${FLOWPACKER_DIR}/config/inference/${kinase}.yaml"
    FINAL_PATH="/store/jaeohshin/work/kine/output/$kinase/final_backbones"

    cp "$BASE_CONFIG" "$KINASE_CONFIG"
    sed -i "s|^  test_path:.*|  test_path: '$FINAL_PATH/'|" "$KINASE_CONFIG"

    pushd "$FLOWPACKER_DIR" > /dev/null
    python sampler_pdb.py "$kinase" "$kinase"
    popd > /dev/null
#    rm "$KINASE_CONFIG"

    echo "✅ Done: $kinase"
done

