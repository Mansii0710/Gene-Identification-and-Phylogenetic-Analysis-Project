#!/usr/bin/python3

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

genomelist = [
    "KP_1705",
    "KP_1899",
    "KP_1904",
    "KP_2146",
    "KP_2343"
]

for genome in genomelist:
    fasta_file = f"{genome}.fasta"
    blast_file = f"{genome}_blast.txt"

    contigs = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    with open(blast_file, "r") as file:
        for line in file:
            columns = line.strip().split("\t")
            query_gene = columns[0]
            contig_id = columns[1]
            perc_identity = float(columns[2])
            sstart = int(columns[8])
            send = int(columns[9])

            if perc_identity >= 95 and contig_id in contigs:
                start = min(sstart, send) - 1 
                end = max(sstart, send)
                sequence = contigs[contig_id].seq[start:end]

                if sstart > send:
                    sequence = sequence.reverse_complement()

                out_file = f"{query_gene}_{genome}.fasta"
                gene_record = SeqRecord(sequence, id=f"{query_gene}_{genome}", description="")
                SeqIO.write(gene_record, out_file, "fasta")
