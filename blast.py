#!/usr/bin/python3

import os

genomelist = [
    "KP_1705",
    "KP_1899",
    "KP_1904",
    "KP_2146",
    "KP_2343"
]

query_file = "Gene_seq.txt"

for genome in genomelist:
    fasta_file = f"{genome}.fasta"
    blast_outfile = f"{genome}_blast.txt"

    print(f"\n Creating database for: {fasta_file}")
    os.system(f"makeblastdb -in {fasta_file} -dbtype nucl")

    print(f"\n Running BLAST for: {fasta_file}")
    os.system(f"blastn -query {query_file} -db {fasta_file} -out {blast_outfile} -outfmt 6")

print("\n All BLAST searches completed.")
