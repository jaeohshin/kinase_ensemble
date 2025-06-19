import os
import glob
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.StructureBuilder import StructureBuilder

def mutate_to_gly_and_renumber(input_pdb, output_pdb):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", input_pdb)

    # Renumbering
    new_resseq = 0
    for model in structure:
        for chain in model:
            new_residues = []
            for residue in list(chain):
                if residue.id[0] == " ":  # Skip waters/heteroatoms
                    residue.id = (" ", new_resseq, " ")
                    residue.resname = "GLY"
                    for atom in list(residue):
                        if atom.get_name() not in ("N", "CA", "C", "O"):
                            residue.detach_child(atom.id)
                    new_resseq += 1

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)

def process_directory(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    pdb_files = sorted(glob.glob(os.path.join(input_dir, "backbone_*.pdb")))

    for pdb_file in pdb_files:
        filename = os.path.basename(pdb_file).replace(".pdb", "_gly.pdb")
        output_file = os.path.join(output_dir, filename)
        print(f"[INFO] Converting {pdb_file} → {output_file}")
        mutate_to_gly_and_renumber(pdb_file, output_file)

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python scripts/mutate_backbone_to_gly.py <bioemu_dir>")
        exit(1)

    input_dir = sys.argv[1]
    output_dir = os.path.join(input_dir, "gly")
    process_directory(input_dir, output_dir)
