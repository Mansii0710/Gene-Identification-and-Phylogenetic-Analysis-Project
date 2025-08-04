# Gene Identification and Phylogenetic Analysis

This project identifies similar genes (gapA, infB, mdh, pgi, phoE, rpoB, and tonB) across multiple genome contigs with:
- 100% query coverage
- Minimum 95% sequence identity

It also builds a phylogenetic tree using Maximum Likelihood with bootstrapping.

## Folder Structure

```
project/
│
├── blast.py                  # Run BLAST between query genes and genome contigs
├── present_genes.py          # Check which genes are present with 100% query coverage and 95% identity 
├── matched_genes.py          # Extract matching gene sequences from genomes
├── paralogs.py               # Check for gene paralogs in each genome
├── conserved_genes.py        # Identify conserved, single-copy genes across genomes
│
├── Gene_seq.txt              # Query gene sequences
├── *.fasta                   # Genome contig files (KP_*.fasta)
├── *_blast.txt               # BLAST output files
├── *_KP_*.fasta              # Extracted gene sequences from genomes
├── *_combined.fasta          # Combined gene sequences for tree building
├── *.meg                     # Aligned files (from MEGA)
├── *.nwk                     # Newick tree files
├── *.png                     # Tree images
```

## Installation 

#Python
sudo apt install python3 python3-pip

#Biopython
pip install biopython

#BLAST
sudo apt install ncbi-blast+

#MEGA
Visit megasoftware.net link and install


## Step-by-Step Execution

### 1. Run BLAST

Then run:

python3 blast.py

This creates databases for each genome and performs BLAST of all 7 genes against them. Output files: *_blast.txt


### 2. Check Which Genes Are Present
Run:

python3 present_genes.py

This script checks for:
- 100% query coverage
- ≥95% identity

It prints which genes are present in each genome. Output: All genes present


### 3. Extract Matching Sequences
Run:

python3 matched_genes.py

This script extracts the matched gene sequences and saves each match in a FASTA file. Output: NO matching genes


### 4. Build Phylogenetic Trees
Steps (manually using MEGA software):
1. Combine all FASTA files of one gene into one file (e.g., `gapA_combined.fasta`) using
   cat gapA_KP_*.fasta > gapA_combined.fasta
2. Align using ClustalW in MEGA and save as `.meg`
3. Build a tree using:
   - Maximum Likelihood method
   - Bootstrap: 100 replicates
4. Export tree as `.nwk` and image `.png`

Repeat for each gene.


### 5. Check for Paralogs
Run:

python3 paralogs.py

This script checks whether a gene appears more than once in each genome. Output: NO paralogs

### 6. Find Other Conserved, Single-Copy Genes
Run:

python3 conserved_genes.py

This identifies genes other than the 7 listed that are:
- Found in all genomes
- Occur only once per genome

Output: NO genes found.



## Answers to Questions

**Q1: Are all 7 genes present in all genomes?**  
Run `present_genes.py` to check presence with 100% coverage and ≥95% identity. Output: All genes present

**Q2: Are there paralogs?**  
Run `paralogs.py`. If a gene appears more than once in any genome, it has paralogs. Output: NO paralogs

**Q3: Which other genes are conserved and unique across genomes?**  
Run `conserved_genes.py`. The output will show how many genes meet this criteria. Output: NO genes found.



## Requirements
- Python 3.6+
- Biopython
- BLAST+ command-line tools (makeblastdb, blastn)
- MEGA software (for alignment and tree building)



