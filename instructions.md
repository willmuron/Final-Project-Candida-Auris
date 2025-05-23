# Downloading data
```
git clone FinalProject
cd FinalProject
```
- You should have a file with the accession numbers called `accession_list.txt'

# Intalling a new version of the sra-tools

```
module load anaconda3
conda create -n sra-tools -c bioconda -c conda-forge sra-tools
```

# Create a script to download the data using the `accession_list.txt` file

```
vi download_sra.sh
```

- Type I
- Copy and past script:

```
#!/bin/bash
# Use prefetch to download all files
module load anaconda3
conda activate sra-tools
while read -r SRR; do
   echo "Downloading $SRR..."
   prefetch --max-size 100G $SRR
   fasterq-dump $SRR --split-files -O fastq_files/
done < accession_list.txt
```
- Type `exit` twice (as needed) to exit the HPC.
- Login again
- go to your ocean folder
- cd into the cloned repo
- Run the download_sra.sh script

```
bash download_sra.sh
```
- This script creates a folder `fastq_files/` which contains all of your sequences.

# Download reference genome for the assembly
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/013/715/GCF_003013715.1_ASM301371v2/GCF_003013715.1_ASM301371v2_genomic.fna.gz
```
```
gunzip GCF_003013715.1_ASM301371v2_genomic.fna.gz
```
```
mv GCF_003013715.1_ASM301371v2_genomic.fna reference.fasta
```

# Build bowtie index
```
module load bowtie2/2.4.4
```
```
bowtie2-build reference.fasta candida_index
```
# Assemble genomes using bowtie2
```
rm assembly.slurm
vi assembly.slurm
```
- Type I, copy and paste the following:
```
#!/bin/bash
#SBATCH --job-name=candida_genome_assembly
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --array=0-14
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/bowtie.out
#SBATCH --mail-user=wmuron@svsu.edu
#SBATCH --mail-type=ALL

#Define samples names
SAMPLES=($(ls fastq_files/*_1.fastq | sed 's|.*/||' | sed 's/_1.fastq//' | sort))
# Get sample name for this task ID
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

echo "Processing $SAMPLE..."

module load bowtie2/2.4.4

bowtie2 --very-fast-local -p 8 -x candida_index -1 fastq_files/${SAMPLE}_1.fastq -2 fastq_files/${SAMPLE}_2.fastq -S ${SAMPLE}.sam

module load samtools
module load bcftools
module load htslib

echo "Converting and Sorting BAM in one step..."
samtools view -@ 8 -bS ${SAMPLE}.sam | samtools sort -@ 8 -m 4G -T /scratch/tmp_sort -o ${SAMPLE}.bam
rm ${SAMPLE}.sam

echo "Indexing BAM file..."
samtools index ${SAMPLE}.bam

echo "Calling SNPs with BCFtools..."
bcftools mpileup -Ou -f reference.fasta ${SAMPLE}.bam | bcftools call -c -Ob -o ${SAMPLE}.bcf

echo "Converting BCF to VCF..."
bcftools view -Ov -o ${SAMPLE}.vcf ${SAMPLE}.bcf

echo "Compressing and indexing VCF..."
bgzip -c ${SAMPLE}.vcf > ${SAMPLE}.vcf.gz
bcftools index ${SAMPLE}.vcf.gz
tabix -p vcf ${SAMPLE}.tab.vcf.gz

echo "Pipeline completed successfully!"
```

```
sbatch assembly.slurm
```
---
# Focusing on antifungal genes:
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/013/715/GCF_003013715.1_ASM301371v2/GCF_003013715.1_ASM301371v2_genomic.gff.gz
gunzip GCF_003013715.1_ASM301371v2_genomic.gff.gz
mv GCF_003013715.1_ASM301371v2_genomic.gff annotation.gff
```
```
grep -i 'CYP51\|FKS1\|1,3-beta-glucan synthase\|CDR1\|Tachykinin' annotation.gff > resistance_genes.gff

awk '$3 == "CDS" {
    match($0, /product=([^;]+)/, a);
    gsub(/%2C/, ",", a[1]);  # decode comma if present
    gsub(/ /, "_", a[1]);   # remove spaces for compatibility
    print $1"\t"$4-1"\t"$5"\t"a[1]
}' resistance_genes.gff > resistance_genes_expanded.bed
```
## Extract antifungal SNPs from vcf files
1. create script to merge vcf files
```
vi merge_vcf.sh
```
-Type I to edit and paste:
```
#!/bin/bash
set -euo pipefail

# Load necessary modules
module load bcftools
module load htslib

# Define directories
VCF_DIR="vcfs"
FILTERED_DIR="filtered_vcfs"
mkdir -p "$FILTERED_DIR"

# Step 1: Extract resistance gene SNPs from each .vcf.gz
for file in "$VCF_DIR"/*.vcf.gz; do
    base=$(basename "$file" .vcf.gz)
    echo "Extracting SNPs from $file..."
    bcftools view -R resistance_genes_expanded.bed "$file" -Oz -o "$FILTERED_DIR/${base}.resistance.vcf.gz"
    tabix -p vcf "$FILTERED_DIR/${base}.resistance.vcf.gz"
done

# Step 2: Merge all filtered VCFs
echo "Merging all resistance VCFs..."
bcftools merge filtered_vcfs/*.resistance.vcf.gz -Oz -o merged_resistance.vcf.gz
tabix -p vcf merged_resistance.vcf.gz
echo "Generating SNP matrix..."
bcftools query -f '%CHROM\t%POS[\t%GT]\n' merged_resistance.vcf.gz > snp_matrix.tsv

echo "All done!"
```
- run script:
```
bash merge_vcf.sh
```

# Extract only informative    SNPS
```
module load bcftools
module load htslib
bcftools view -v snps merged_resistance.vcf.gz -Oz -o snps_only.vcf.gz
tabix -p vcf snps_only.vcf.gz
bedtools intersect -a snps_only.vcf.gz -b resistance_genes_expanded.bed -wa -wb > snps_with_genes.txt
bcftools query -l snps_only.vcf.gz > sample_names.txt
```
- push to repo

# Visualize snps in antifungal genes in R
```R
library(tidyverse)

# Load SNP + gene file
df <- read.table("snps_with_genes2.txt", header = FALSE, stringsAsFactors = FALSE)

# Load sample names and metadata
samples <- readLines("sample_names.txt")
metadata <- read.table("sample_metadata.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Assign column names
colnames(df)[1:9] <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
colnames(df)[10:(9 + length(samples))] <- samples
colnames(df)[(10 + length(samples)):(10 + length(samples) + 2)] <- c("GENE_CHROM", "GENE_START", "GENE_END")
colnames(df)[(10 + length(samples) + 3)] <- "GENE"

# Extract relevant columns and reshape
genos <- df[, c("CHROM", "POS", samples, "GENE")]

genos_long <- genos %>%
  pivot_longer(cols = all_of(samples), names_to = "Sample", values_to = "GT_RAW") %>%
  mutate(GT = sub(":.*", "", GT_RAW)) %>%
  filter(GT %in% c("0/1", "1/1")) %>%
  count(GENE, Sample)

metadata <- metadata %>%
  mutate(
    SampleID = sub("\\.bam$", "", SamplsID),  # strip ".bam"
    Label = paste(SampleID, Country, Source, sep = "_")
  )

annotated <- genos_long %>%
  mutate(Sample = sub("\\.bam$", "", Sample)) %>%  # strip ".bam" from genos_long
  left_join(metadata, by = c("Sample" = "SampleID"))

library(viridis)

ggplot(annotated, aes(x = Label, y = GENE, fill = n)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradientn(
    colors = c("navy", "white", "firebrick"),
    values = scales::rescale(c(min(annotated$n), median(annotated$n), max(annotated$n))),
    name = "SNP count"
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "SNPs per Gene per Sample",
    x = NULL,
    y = "Gene"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.grid = element_blank()
  )
```

# PHYLOGENETIC TREE

# Step 3: Convert merged VCF to SNP matrix

# Let convert vcf to phyllip to run phylogenetic analysis (avoids mafft) - Make sure vcf2phylip file is in the current directory
```
module load python/3.8.6 
python vcf2phylip.py -i merged_resistance.vcf.gz -m2
```
- NOTE that the python script vcf2phylip.py has been added to this repository and will be downloaded to the HPC when you clone it.

- Clean file to remove samples with no nucleotides

```
awk 'NR==1 {len=$2; next} 
     $2 ~ /[^Nn?-]/ {print $0 > "filtered_body.phy"; c++} 
     END {print c, len > "filtered_head.phy"}' merged_resistance.min2.phy

cat filtered_head.phy filtered_body.phy > cleaned.phy
rm filtered_head.phy filtered_body.phy
```

# Run Raxml
```
module load RAxML
raxmlHPC -s cleaned.phy -n antifungal3 -m GTRCAT -p 12345 -x 12345 -# 1000 -f a
```
- change merged_resistance.min2.phy with name of file produced by python script above

# Visualize tree in iTOL
https://itol.embl.de/upload.cgi

- Make sure you display bootstrap support value in nodes.




