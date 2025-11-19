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
NUM_FRAMES="${2:-50}"

INPUT_DIR="/store/jaeohshin/tools/flowpacker/samples/${PROTEIN_NAME}/run_1"
OUTPUT_DIR="/store/jaeohshin/work/kine/output/${PROTEIN_NAME}/final_str"
WORK_DIR="/store/jaeohshin/work/kine/output/${PROTEIN_NAME}/gromacs_work"

EM_MDP="/store/jaeohshin/work/kine/scripts/em.mdp"
NVT_MDP="/store/jaeohshin/work/kine/scripts/nvt.mdp"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

process_frame() {
    local i=$(printf "%03d" $1)
    echo "=== Processing frame $i ==="

    INPUT_PDB="${INPUT_DIR}/final_${i}.pdb"
    OUTPUT_PDB="${OUTPUT_DIR}/receptor_${i}.pdb"
    if [ -f "$OUTPUT_PDB" ]; then echo "✓ Already processed"; return 0; fi
    if [ ! -f "$INPUT_PDB" ]; then echo "✗ Missing input"; return 1; fi

    echo "Step 1: Generate processed PDB + topology"
    gmx pdb2gmx -f "$INPUT_PDB" -ignh -o prot_${i}.pdb -p topol_${i}.top <<EOF
14
1
EOF

    echo "→ Creating heavy atom restraints (EM)"
    gmx select -f prot_${i}.pdb -s prot_${i}.pdb -select 'not name "H*"' -on heavy_${i}.ndx
    gmx genrestr -f prot_${i}.pdb -n heavy_${i}.ndx -o posre_heavy_${i}.itp -fc 500 500 500

    echo "→ Creating Cα restraints (NVT)"
    gmx select -f prot_${i}.pdb -s prot_${i}.pdb -select 'name CA' -on ca_${i}.ndx
    gmx genrestr -f prot_${i}.pdb -n ca_${i}.ndx -o posre_Calpha_${i}.itp -fc 5 5 5

    # Replace default posre.itp with heavy restraints for EM
    sed -i "s/posre.itp/posre_heavy_${i}.itp/" topol_${i}.top

    echo "Step 2: Box + solvate + ionize"
    gmx editconf -f prot_${i}.pdb -o prot_box_${i}.gro -d 1.0 -bt cubic -c

#    gmx solvate -cp prot_box_${i}.pdb -cs -o prot_solv_${i}.pdb -p topol_${i}.top
#    gmx grompp -f "$EM_MDP" -c prot_solv_${i}.pdb -p topol_${i}.top -r prot_solv_${i}.pdb -o genion_${i}.tpr -maxwarn 40
    gmx solvate -cp prot_box_${i}.gro -cs -o prot_solv_${i}.gro -p topol_${i}.top
    gmx grompp -f "$EM_MDP" -c prot_solv_${i}.gro -p topol_${i}.top -r prot_solv_${i}.gro -o genion_${i}.tpr -maxwarn 40

    gmx genion -s genion_${i}.tpr -neutral -conc 0.15 -o prot_solv_ions_${i}.gro -p topol_${i}.top <<EOF
SOL
EOF

    echo "Step 3: EM"
    gmx grompp -f "$EM_MDP" -c prot_solv_ions_${i}.gro -p topol_${i}.top -r prot_solv_ions_${i}.gro -o minim_${i}.tpr -maxwarn 40
    gmx mdrun -s minim_${i}.tpr -nt 8  -deffnm minim_${i}

    gmx trjconv -s minim_${i}.tpr -f minim_${i}.gro -pbc mol -ur compact -center -o minim_whole_${i}.gro <<EOF
Protein
System
EOF

    echo "Step 4: Switch to Cα restraints for NVT"
    sed -i "s/posre_heavy_${i}.itp/posre_Calpha_${i}.itp/" topol_${i}.top
    gmx grompp -f "$NVT_MDP" -c minim_whole_${i}.gro -p topol_${i}.top -r minim_whole_${i}.gro -o nvt_${i}.tpr -maxwarn 40
    gmx mdrun -s nvt_${i}.tpr -nt 8 -deffnm nvt_${i}

    gmx make_ndx -f nvt_${i}.gro -o protein_${i}.ndx <<EOF
1
q
EOF
    gmx editconf -f nvt_${i}.gro -o temp_final_${i}.pdb -n protein_${i}.ndx <<EOF
1
EOF
    cp temp_final_${i}.pdb "$OUTPUT_PDB"
    echo "✓ Saved → $OUTPUT_PDB"

    # Optional cleanup
    rm -f heavy_${i}.ndx ca_${i}.ndx posre_heavy_${i}.itp posre_Calpha_${i}.itp
}


for i in $(seq 0 $((NUM_FRAMES-1))); do
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

