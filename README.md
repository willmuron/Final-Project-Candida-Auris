# Final-Project
# Phylogenetic Analysis of Antifungal Resistance Genes in *Candida auris*

## 🧬 Project Overview
This project investigates the evolutionary relationships among antifungal resistance genes in the fungal pathogen *Candida auris*. We focus on genes such as **ERG11**, **FKS1**, and **CDR1**, which are commonly associated with resistance to antifungal medications. Using bioinformatics tools, we perform multiple sequence alignment and construct phylogenetic trees to identify mutations and clustering patterns.

## 🧠 Biological Problem
*Candida auris* is an emerging, multidrug-resistant fungal pathogen responsible for hospital outbreaks and invasive infections. The rapid global spread and difficulty in treatment make it a public health concern. By studying the evolution of its resistance genes, we aim to better understand its resistance mechanisms and potential transmission pathways.

"How have antifungal resistance genes in Candida auris evolved, and what mutations are associated with resistance in different clinical isolates?"

## 🧪 Hypothesis
We hypothesize that resistant strains of *C. auris* will:
- Contain shared mutations in key resistance genes.
- Cluster together in phylogenetic trees based on sequence similarity.

## 🎯 Goals and Objectives
- Collect antifungal resistance gene sequences from publicly available *C. auris* genomes.
- Perform multiple sequence alignment (MSA) using MAFFT or MUSCLE.
- Construct phylogenetic trees using RAxML or FastTree.
- Identify clustering patterns and evolutionary trends associated with resistance.

## 🛠️ Tools and Technologies
- [MAFFT](https://mafft.cbrc.jp/alignment/software/) – Sequence alignment
- [MUSCLE](https://www.drive5.com/muscle/) – Sequence alignment
- [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) – Gene identification
- [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/) – Phylogenetic tree construction
- [FastTree](http://www.microbesonline.org/fasttree/) – Alternative tree-building method
- [FigTree](https://github.com/rambaut/figtree) / [iTOL](https://itol.embl.de/) – Tree visualization
- Python / R – Data processing and plotting

- ### 🧰 HPC Software Availability
We verified the following tools are available on the MSU HPC:

- MAFFT
- MUSCLE
- BLAST+
- FastQC
- FastTree
- sratoolkit

These tools will be used for alignment, gene identification, quality control, and tree construction.

### 🔧 Tools to Be Installed
We may need to install:
- FigTree (used locally for tree visualization)
- iTOL (web-based, no install required)


## 📋 How to Run the Analysis

### 📂 SRA Datasets
We are using publicly available *Candida auris* genome data from the NCBI Sequence Read Archive (SRA). Accession numbers are listed in `accession_list.txt`.

We selected 15 isolates from different geographic regions and clades, including samples from the United States, United Kingdom, Colombia, Iran, and South Africa. The list also includes representatives from key phylogenetic clades (South Asian, South American, and South African) and samples that are non-clinical in nature.

To download the datasets, use:
```bash
prefetch --option-file accession_list.txt

## 🔬 Project Workflow

This section outlines the step-by-step procedure we are using to complete our analysis of antifungal resistance genes in *Candida auris*.

---

### ✅ Current Status
- ✔️ Project topic defined
- ✔️ Data collection plan created
- ✔️ `accession_list.txt` with SRR IDs completed
- ✔️ Project README structured and documented

---

### 🔄 Next Steps

#### **1. Download Raw Data**
Use the NCBI SRA Toolkit to download the raw FASTQ data:

```bash
# Download using accession list
prefetch --option-file accession_list.txt

# Convert .sra files to .fastq
for sra in $(cat accession_list.txt); do
    fastq-dump --split-files $sra
done

# ===============================
# STEP 1: Download Raw FASTQ Data
# ===============================

# Download all SRA files listed in accession_list.txt
prefetch --option-file accession_list.txt

# Convert each .sra file to paired-end FASTQ files
for sra in $(cat accession_list.txt); do
    fastq-dump --split-files $sra
done

# ===============================
# STEP 2: Run Quality Control
# ===============================

# Create a directory to hold FastQC reports
mkdir -p fastqc_reports

# Run FastQC on all downloaded FASTQ files
fastqc *.fastq -o fastqc_reports

# ===============================
# STEP 3: (Optional) Trim Reads
# ===============================

# Example Trimmomatic command for paired-end trimming
# (Run only if FastQC reports indicate issues)
trimmomatic PE SRRXXXX_1.fastq SRRXXXX_2.fastq \
  SRRXXXX_1_trimmed.fastq SRRXXXX_1_unpaired.fastq \
  SRRXXXX_2_trimmed.fastq SRRXXXX_2_unpaired.fastq \
  SLIDINGWINDOW:4:20 MINLEN:50

# ===============================
# STEP 4: Extract Resistance Genes
# ===============================

# OPTION A: If working with full genomes (BLAST approach)
makeblastdb -in genome.fasta -dbtype nucl
blastn -query ERG11_reference.fasta -db genome.fasta -out ERG11_hits.txt

# OPTION B: If working with raw reads (mapping + extraction approach)

# Index the reference genome
bwa index reference_genome.fasta

# Map reads to reference genome
bwa mem reference_genome.fasta SRRXXXX_1.fastq SRRXXXX_2.fastq > mapped_reads.sam

# Convert to BAM and sort
samtools view -Sb mapped_reads.sam | samtools sort -o sorted_reads.bam

# Index the sorted BAM
samtools index sorted_reads.bam

# Extract the gene region (replace "chr:start-end" with actual coordinates)
samtools faidx reference_genome.fasta chr:start-end > ERG11_sequence.fasta

# ===============================
# STEP 5: Multiple Sequence Alignment
# ===============================

# Align gene sequences using MAFFT
mafft --auto resistance_genes.fasta > aligned_genes.fasta

# ===============================
# STEP 6: Phylogenetic Tree Construction
# ===============================

# OPTION A: Use FastTree (quick method)
FastTree -nt aligned_genes.fasta > resistance_tree.nwk

# OPTION B: Use RAxML (maximum likelihood + bootstrap)
raxmlHPC -s aligned_genes.fasta -n resistance_tree -m GTRGAMMA -p 12345 -# 100

# ===============================
# OUTPUT: Tree file is in Newick format (.nwk)
# Visualize using iTOL or FigTree
# ===============================

