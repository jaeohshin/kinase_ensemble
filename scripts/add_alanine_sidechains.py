import os
import glob
from Bio.PDB import PDBParser, PDBIO


def mutate_to_ala_and_renumber(input_pdb, output_pdb):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", input_pdb)

    new_resseq = 1
    for model in structure:
        for chain in model:
            new_residues = []
            for residue in list(chain):
                if residue.id[0] != " ":  # Skip hetero/water
                    continue

                # Make a copy of residue with filtered atoms
                atoms_to_keep = []
                for atom in residue:
                    if atom.get_name() in ("N", "CA", "CB", "C", "O"):
                        atoms_to_keep.append(atom.copy())

                # Create new residue
                from Bio.PDB import Residue
                new_residue = Residue.Residue((" ", new_resseq, " "), "ALA", residue.segid)
                for atom in atoms_to_keep:
                    new_residue.add(atom)
                new_residues.append(new_residue)
                new_resseq += 1

            # Remove all residues from chain
            for residue in list(chain):
                chain.detach_child(residue.id)

            # Add new clean residues
            for residue in new_residues:
                chain.add(residue)

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)


def process_directory(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    pdb_files = sorted(glob.glob(os.path.join(input_dir, "backbone_*.pdb")))

    for pdb_file in pdb_files:
        filename = os.path.basename(pdb_file).replace(".pdb", "_ala.pdb")
        output_file = os.path.join(output_dir, filename)
        print(f"[INFO] Converting {pdb_file} → {output_file}")
        mutate_to_ala_and_renumber(pdb_file, output_file)

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python scripts/mutate_backbone_to_ala.py <bioemu_dir>")
        exit(1)

    input_dir = sys.argv[1]
    output_dir = os.path.join(input_dir, "ala")
    process_directory(input_dir, output_dir)

