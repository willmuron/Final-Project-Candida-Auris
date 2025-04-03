# Phylogenetic Analysis of Antifungal Resistance Genes in *Candida auris*

## 🧬 Project Overview
This project investigates the evolutionary relationships among antifungal resistance genes in the fungal pathogen *Candida auris*. We focus on genes such as **ERG11**, **FKS1**, and **CDR1**, which are commonly associated with resistance to antifungal medications. Using bioinformatics tools, we perform multiple sequence alignment and construct phylogenetic trees to identify mutations and clustering patterns.

## 🧠 Biological Problem
*Candida auris* is an emerging, multidrug-resistant fungal pathogen responsible for hospital outbreaks and invasive infections. The rapid global spread and difficulty in treatment make it a public health concern. By studying the evolution of its resistance genes, we aim to better understand its resistance mechanisms and potential transmission pathways.

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

---

## 🔬 Analysis Workflow

### ✅ Step 1: Define Biological Question
This step is outlined above in the **Biological Problem** and **Hypothesis** sections.

---

### 📥 Step 2: Download Raw Data

```bash
# Download all SRA files listed in accession_list.txt
prefetch --option-file accession_list.txt

# Convert each .sra file to paired-end FASTQ files
for sra in $(cat accession_list.txt); do
    fastq-dump --split-files $sra
done
```

---

### ✅ Step 3: Quality Control

```bash
# Create a directory to hold FastQC reports
mkdir -p fastqc_reports

# Run FastQC on all downloaded FASTQ files
fastqc *.fastq -o fastqc_reports
```

---

### ✂️ Step 4: (Optional) Read Trimming

```bash
# Example Trimmomatic command for paired-end trimming
# (Edit and uncomment if needed)
# trimmomatic PE SRRXXXX_1.fastq SRRXXXX_2.fastq \
#   SRRXXXX_1_trimmed.fastq SRRXXXX_1_unpaired.fastq \
#   SRRXXXX_2_trimmed.fastq SRRXXXX_2_unpaired.fastq \
#   SLIDINGWINDOW:4:20 MINLEN:50
```

---

### 🧬 Step 5: Gene Extraction or Mapping

**Option A – Full genomes with BLAST:**

```bash
makeblastdb -in genome.fasta -dbtype nucl
blastn -query ERG11_reference.fasta -db genome.fasta -out ERG11_hits.txt
```

**Option B – Raw reads mapped to reference:**

```bash
# Index the reference genome
bwa index reference_genome.fasta

# Map reads
bwa mem reference_genome.fasta SRRXXXX_1.fastq SRRXXXX_2.fastq > mapped_reads.sam

# Convert to sorted BAM
samtools view -Sb mapped_reads.sam | samtools sort -o sorted_reads.bam
samtools index sorted_reads.bam

# Extract gene region (replace chr:start-end with actual coordinates)
samtools faidx reference_genome.fasta chr:start-end > ERG11_sequence.fasta
```

---

### 🧬 Step 6: Multiple Sequence Alignment

```bash
mafft --auto resistance_genes.fasta > aligned_genes.fasta
```

---

### 🌲 Step 7: Phylogenetic Tree Construction

**Option A – FastTree:**

```bash
FastTree -nt aligned_genes.fasta > resistance_tree.nwk
```

**Option B – RAxML:**

```bash
raxmlHPC -s aligned_genes.fasta -n resistance_tree -m GTRGAMMA -p 12345 -# 100
```

---

### 🌳 Step 8: Tree Visualization

Use one of the following tools:
- [iTOL (web-based)](https://itol.embl.de/)
- [FigTree (desktop-based)](https://github.com/rambaut/figtree)

Upload the `resistance_tree.nwk` file to view and annotate your tree.

---

### 🧬 Step 9: (Optional) Mutation Analysis

Use tools like:
- [MEGA](https://www.megasoftware.net/)
- [Geneious](https://www.geneious.com/)
- [SnpEff](https://pcingola.github.io/SnpEff/)

To inspect and annotate specific mutations in aligned sequences.

---

## ✅ Final Output
- Aligned resistance gene sequences: `aligned_genes.fasta`
- Phylogenetic tree: `resistance_tree.nwk`
- FastQC reports
- Optional: mutation annotations and figures
