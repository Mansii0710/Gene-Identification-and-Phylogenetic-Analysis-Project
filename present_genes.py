#!/usr/bin/python3

from Bio import SeqIO

genomelist = [
    "KP_1705.fasta",
    "KP_1899.fasta",
    "KP_1904.fasta",
    "KP_2146.fasta",
    "KP_2343.fasta"
]

gene_len = {}
for gene in SeqIO.parse("Gene_seq.txt", "fasta"):
    gene_len[gene.id] = len(gene.seq)

for genome in genomelist:
    out_file = genome.replace(".fasta", "_blast.txt")
    matched = []

    print(f"\n {out_file}")

    with open(out_file, "r") as file:
        for line in file:
            columns = line.strip().split("\t")
            query_gene = columns[0]
            perc_identity = float(columns[2])
            align_len = int(columns[3])

            if query_gene in gene_len:
                coverage = (align_len / gene_len[query_gene]) * 100
                if align_len == gene_len[query_gene] and perc_identity >= 95:
                    print(f"{query_gene}   Query Coverage: {coverage:.2f}%, Identity: {perc_identity:.2f}%")
                    if query_gene not in matched:
                        matched.append(query_gene)

    for gene in gene_len:
        if gene in matched:
            print(f"{gene}: Present")
        else:
            print(f"{gene}: Not found")
