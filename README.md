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

## 📁 Project Structure
