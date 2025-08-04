#!/usr/bin/python3

from collections import defaultdict
from Bio import SeqIO

genomelist = [
    "KP_1705",
    "KP_1899",
    "KP_1904",
    "KP_2146",
    "KP_2343"
]

check_genes = {"gapA", "infB", "mdh", "pgi", "phoE", "rpoB", "tonB"}

gene_len = {gene.id: len(gene.seq) for gene in SeqIO.parse("Gene_seq.txt", "fasta")}


for genome in genomelist:
    blast_file = f"{genome}_blast.txt"
    gene_hits = defaultdict(int)

    with open(blast_file, "r") as file:
        for line in file:
            cols = line.strip().split("\t")
            query_gene = cols[0]
            perc_identity = float(cols[2])
            align_len = int(cols[3])

            if query_gene in check_genes:
                coverage = (align_len / gene_len[query_gene]) * 100
                if perc_identity >= 95 and coverage >= 100:
                    gene_hits[query_gene] += 1

    for gene, count in gene_hits.items():
        if count > 1:
            print(f"{gene} appears {count} times in {genome}")
