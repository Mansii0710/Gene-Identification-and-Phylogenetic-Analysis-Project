#!/usr/bin/python3

from collections import defaultdict
from Bio import SeqIO

genomelist = ["KP_1705", "KP_1899", "KP_1904", "KP_2146", "KP_2343"]
unwanted_geneslist = {"gapA", "infB", "mdh", "pgi", "phoE", "rpoB", "tonB"}

gene_len = {gene.id: len(gene.seq) for gene in SeqIO.parse("Gene_seq.txt", "fasta")}

hit_counts = defaultdict(lambda: defaultdict(int)) 

for genome in genomelist:
    with open(f"{genome}_blast.txt", "r") as file:
        for line in file:
            parts = line.strip().split("\t")
            query_gene = parts[0]
            perc_identity = float(parts[2])
            align_len = int(parts[3])

            if query_gene in unwanted_geneslist:
                continue

            if query_gene in gene_len:
                coverage = (align_len / gene_len[query_gene]) * 100
                if perc_identity >= 95 and coverage == 100:
                    hit_counts[query_gene][genome] += 1

conserved_genes = []

for gene, genome_hits in hit_counts.items():
    if all(genome in genome_hits and genome_hits[genome] == 1 for genome in genomelist):
        conserved_genes.append(gene)

print("\n Conserved and unique genes across all genomes:")
for gene in conserved_genes:
    print(f"- {gene}")

print(f"\n Total genes found: {len(conserved_genes)}")
