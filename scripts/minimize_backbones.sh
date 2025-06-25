#!/bin/bash

# Usage: bash scripts/strip_and_minimize.sh output/abl1/bioemu/gly

INPUT_DIR=$1
WATERMODEL="none"         # vacuum
FORCEFIELD_ID=14          # AMBER99SB-ILDN

if [ -z "$INPUT_DIR" ]; then
    echo "Usage: bash scripts/strip_and_minimize.sh <glycine_structs_dir>"
    exit 1
fi

cd "$INPUT_DIR" || exit

for pdb in backbone_*_gly.pdb; do
    id="${pdb%.pdb}"
    echo "[INFO] Processing $pdb"

    # Skip if minimized file already exists
    if [ -f "minimized_${pdb}" ]; then
        echo "[INFO] minimized_${pdb} already exists, skipping."
        continue
    fi

    # Step 1: Clean previous files for this structure
    rm -f em_${id}.* ${id}_processed.gro ${id}_boxed.gro topol.top \#* *.itp *.log *.edr *.trr

    # Step 2: Generate topology and add hydrogens
    gmx pdb2gmx -f "$pdb" -o "${id}_processed.gro" -p topol.top -water "$WATERMODEL" <<EOF
$FORCEFIELD_ID
1
1
EOF

    # Lower restraint constants in posre.itp (optional step)
    # Change only the last three fields (fx fy fz) to 50
    awk '
    /^[[:space:]]*[0-9]+[[:space:]]+1[[:space:]]+[0-9.]+[[:space:]]+[0-9.]+[[:space:]]+[0-9.]+/ {
    $3 = $4 = $5 = 50
    }
    { print }
    ' posre.itp > posre_fixed.itp && mv posre_fixed.itp posre.itp



    if [ ! -f "${id}_processed.gro" ]; then
        echo "[ERROR] Failed to generate ${id}_processed.gro — skipping."
        continue
    fi

    # Step 3: Create box (vacuum, so just cubic space)
    gmx editconf -f "${id}_processed.gro" -o "${id}_boxed.gro" -c -d 1.0 -bt cubic

    if [ ! -f "${id}_boxed.gro" ]; then
        echo "[ERROR] Failed to create boxed structure — skipping."
        continue
    fi

    # Step 4: Create MDP
cat > em_${id}.mdp <<EOF

define        = -DPOSRES
integrator    = steep
emtol         = 100.0
emstep        = 0.01
nsteps        = 5000

cutoff-scheme = Verlet
coulombtype   = cutoff
rcoulomb      = 1.0
vdwtype       = cutoff
rvdw          = 1.0

constraints   = none
pbc           = xyz
EOF


    # Step 5: grompp
    gmx grompp -f em_${id}.mdp -c "${id}_boxed.gro" -r "${id}_boxed.gro" -p topol.top -o em_${id}.tpr
    if [ $? -ne 0 ]; then
        echo "[ERROR] grompp failed for $id — skipping."
        continue
    fi

    # Step 6: Run minimization
    gmx mdrun -v -deffnm em_${id}

    # Step 7: Save minimized structure
    gmx trjconv -s em_${id}.tpr -f em_${id}.gro -o "minimized_${pdb}" <<< 0
    echo "[INFO] Saved: minimized_${pdb}"

done
