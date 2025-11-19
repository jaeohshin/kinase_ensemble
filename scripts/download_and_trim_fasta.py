#!/usr/bin/env python3

import sys
import os
import requests
from Bio import SeqIO
from io import StringIO

def download_pdb(pdbid):
    url = f"https://files.rcsb.org/download/{pdbid}.pdb"
    response = requests.get(url)
    if not response.ok:
        raise Exception(f"Failed to download PDB: {pdbid}")
    return response.text

def extract_seqres_fasta(pdb_text):
    pdb_io = StringIO(pdb_text)
    return next(SeqIO.parse(pdb_io, "pdb-seqres"))

def main(tsv_file):
    with open(tsv_file) as f:
        for line in f:
            if not line.strip():
                continue
            kinase, pdbid = line.strip().split()
            kinase_lower = kinase.lower()
            print(f"[INFO] Downloading SEQRES for {kinase} ({pdbid})")

            try:
                pdb_text = download_pdb(pdbid)
                fasta_record = extract_seqres_fasta(pdb_text)
                fasta_record.id = kinase_lower
                fasta_record.description = f"from PDB {pdbid}"

                out_name = f"{kinase_lower}.fasta"
                with open(out_name, "w") as out_f:
                    SeqIO.write(fasta_record, out_f, "fasta")
                print(f"[OK] Saved: {out_name}")

            except Exception as e:
                print(f"[ERROR] {kinase}: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python download_fasta_seqres_only.py kinase.txt")
        sys.exit(1)
    main(sys.argv[1])

