#!/bin/bash
set -e

usage() {
    echo "Usage: $0 <protein_name> [num_frames]"
    echo "Example: $0 abl1 100"
    exit 1
}

if [ -z "$1" ]; then
    echo "ERROR: Protein name not specified"
    usage
fi

PROTEIN_NAME="$1"
NUM_FRAMES="${2:-100}"

INPUT_DIR="/store/jaeohshin/tools/flowpacker/samples/${PROTEIN_NAME}/run_1"
OUTPUT_DIR="/store/jaeohshin/work/kine/output/${PROTEIN_NAME}/final_str"
WORK_DIR="/store/jaeohshin/work/kine/output/${PROTEIN_NAME}/gromacs_work"

EM_MDP="/store/jaeohshin/work/kine/scripts/em.mdp"
NVT_MDP="/store/jaeohshin/work/kine/scripts/nvt.mdp"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

# One-time topology and restraints generation
prepare_topology_once() {
    if [ -f "topol.top" ]; then
        echo "✓ Reusing existing topol.top"
        return
    fi

    FIRST_FRAME="${INPUT_DIR}/final_001.pdb"
    if [ ! -f "$FIRST_FRAME" ]; then
        echo "ERROR: $FIRST_FRAME not found"
        exit 1
    fi

    echo "→ Generating topology and restraints from $FIRST_FRAME"
    gmx pdb2gmx -f "$FIRST_FRAME" -ignh -o ref_prot.pdb -p topol.top <<EOF
14
1
EOF

    echo "→ Creating heavy atom restraints (EM)"
    gmx select -f ref_prot.pdb -s ref_prot.pdb -select 'not name "H*"' -on heavy.ndx
    gmx genrestr -f ref_prot.pdb -n heavy.ndx -o posre_heavy.itp -fc 1000 1000 1000

    echo "→ Creating Cα restraints (NVT)"
    gmx select -f ref_prot.pdb -s ref_prot.pdb -select 'name CA' -on ca.ndx
    gmx genrestr -f ref_prot.pdb -n ca.ndx -o posre_Calpha.itp -fc 10 10 10

    rm -f heavy.ndx ca.ndx ref_prot.pdb
}

# Per-frame processing
process_frame() {
    local i=$(printf "%03d" $1)
    echo "=== Processing frame $i ==="

    INPUT_PDB="${INPUT_DIR}/final_${i}.pdb"
    OUTPUT_PDB="${OUTPUT_DIR}/receptor_${i}.pdb"
    if [ -f "$OUTPUT_PDB" ]; then echo "✓ Already processed"; return 0; fi
    if [ ! -f "$INPUT_PDB" ]; then echo "✗ Missing input"; return 1; fi

    cp topol.top topol_${i}.top

    echo "Step 1: Generate processed PDB only"
    gmx pdb2gmx -f "$INPUT_PDB" -ignh -o prot_${i}.pdb -p dummy.top <<EOF
14
1
EOF
    rm -f dummy.top

    echo "Step 2: Set EM restraints"
    sed -i "s/posre.itp/posre_heavy.itp/" topol_${i}.top

    gmx editconf -f prot_${i}.pdb -o prot_box_${i}.pdb -d 1.0 -bt cubic -c
    gmx solvate -cp prot_box_${i}.pdb -cs -o prot_solv_${i}.pdb -p topol_${i}.top
    gmx grompp -f "$EM_MDP" -c prot_solv_${i}.pdb -p topol_${i}.top -r prot_solv_${i}.pdb -o genion_${i}.tpr -maxwarn 40
    gmx genion -s genion_${i}.tpr -neutral -conc 0.15 -o prot_solv_ions_${i}.pdb -p topol_${i}.top <<EOF
SOL
EOF

    gmx grompp -f "$EM_MDP" -c prot_solv_ions_${i}.pdb -p topol_${i}.top -r prot_solv_ions_${i}.pdb -o minim_${i}.tpr -maxwarn 40
    gmx mdrun -s minim_${i}.tpr -gpu_id 0 -nt 16 -deffnm minim_${i}

    gmx trjconv -s minim_${i}.tpr -f minim_${i}.gro -pbc mol -ur compact -center -o minim_whole_${i}.pdb <<EOF
Protein
System
EOF

    echo "Step 3: Switch to Cα restraints"
    sed -i "s/posre_heavy.itp/posre_Calpha.itp/" topol_${i}.top
    gmx grompp -f "$NVT_MDP" -c minim_whole_${i}.pdb -p topol_${i}.top -r prot_solv_ions_${i}.pdb -o nvt_${i}.tpr -maxwarn 40
    gmx mdrun -s nvt_${i}.tpr -gpu_id 0 -nt 16 -deffnm nvt_${i}

    gmx make_ndx -f nvt_${i}.gro -o protein_${i}.ndx <<EOF
a Protein
q
EOF
    gmx editconf -f nvt_${i}.gro -o temp_final_${i}.pdb -n protein_${i}.ndx <<EOF
Protein
EOF
    cp temp_final_${i}.pdb "$OUTPUT_PDB"
    echo "✓ Saved → $OUTPUT_PDB"

    #rm -f *#* prot_${i}.pdb prot_box_${i}.pdb prot_solv_${i}.pdb prot_solv_ions_${i}.pdb minim_whole_${i}.pdb
    #rm -f genion_${i}.tpr minim_${i}.tpr nvt_${i}.tpr
    #rm -f topol_${i}.top protein_${i}.ndx temp_final_${i}.pdb
}

prepare_topology_once

for i in $(seq 1 $NUM_FRAMES); do
    if process_frame $i; then
        echo "Frame $i: SUCCESS"
    else
        echo "Frame $i: FAILED"
    fi
    echo ""
done

success_count=$(ls ${OUTPUT_DIR}/receptor_*.pdb 2>/dev/null | wc -l)
echo "✓ Processed: ${success_count}/${NUM_FRAMES} frames"
echo "✓ Output: $OUTPUT_DIR"
