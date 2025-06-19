from Bio import SeqIO
from Bio.Data.IUPACData import protein_letters_1to3
import os
import glob
from collections import defaultdict

# === Paths ===
fasta_path = "/store/jaeohshin/work/kine/output/abl1/bioemu/sequence.fasta"
input_dir = "/store/jaeohshin/work/kine/output/abl1/bioemu/gly"
output_dir = "/store/jaeohshin/work/kine/output/abl1/bioemu/final_backbones"

# === Load sequence ===
seq_record = next(SeqIO.parse(fasta_path, "fasta"))
aa_seq = str(seq_record.seq)
resnames = [protein_letters_1to3[aa].upper() for aa in aa_seq]

# === Create output dir ===
os.makedirs(output_dir, exist_ok=True)

# === Process each input PDB ===
pdb_files = sorted(glob.glob(os.path.join(input_dir, "minimized_backbone_???_gly.pdb")))

for pdb_file in pdb_files:
    corrected_lines = []
    residue_blocks = defaultdict(list)

    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM"):
                chain_id = line[21]
                res_id = int(line[22:26])
                res_uid = (chain_id, res_id)
                residue_blocks[res_uid].append(line)
            else:
                corrected_lines.append(line)  # preserve headers like REMARK/CRYST1/etc.

    # Replace residue names and renumber from 1
    new_res_idx = 1
    for i, (res_uid, atom_lines) in enumerate(sorted(residue_blocks.items(), key=lambda x: x[0][1])):
        if i >= len(resnames):
            raise ValueError(f"Too many residues in {pdb_file}")
        new_resname = resnames[i]
        for line in atom_lines:
            new_line = (
                line[:17] + f"{new_resname:>3}" +   # residue name
                line[20:22] + f"{new_res_idx:>4}" + # new residue ID
                line[26:]
            )
            corrected_lines.append(new_line)
        new_res_idx += 1

    if new_res_idx - 1 != len(resnames):
        raise ValueError(f"Residue count mismatch in {pdb_file}: got {new_res_idx-1}, expected {len(resnames)}")

    # Write to clean output dir
    base = os.path.basename(pdb_file).replace("_gly.pdb", "_corrected.pdb")
    out_path = os.path.join(output_dir, base)

    with open(out_path, "w") as out_f:
        for line in corrected_lines:
            if not (line.startswith("TER") or line.startswith("ENDMDL")):
                out_f.write(line)
        out_f.write("TER\n")
        out_f.write("ENDMDL\n")
        

print(f"✅ Corrected {len(pdb_files)} files to: {output_dir}")
