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
vi assembly.slurm
```
- Type I, copy and paste the following:
```
#!/bin/bash
#SBATCH --job-name=candida_genome_assembly
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
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

bowtie2 --very-fast-local -p 16 -x candida_index -1 fastq_files/${SAMPLE}_1.fastq -2 fastq_files/${SAMPLE}_2.fastq -S ${SAMPLE}.sam

module load samtools

echo "Converting and Sorting BAM in one step..."
samtools view -@ 16 -bS ${SAMPLE}.sam | samtools sort -@ 16 -m 16G -T /scratch/tmp_sort -o ${SAMPLE}.bam
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
grep -i 'ERG11\|FKS1\|FKS2\|CDR1\|TAC1' annotation.gff > resistance_genes.gff

awk '$3 == "CDS" {
    match($0, /product=([^;]+)/, a);
    gsub(/%2C/, ",", a[1]);  # decode comma if present
    gsub(/ /, "_", a[1]);   # remove spaces for compatibility
    print $1"\t"$4-1"\t"$5"\t"a[1]
}' resistance_genes.gff > resistance_genes.bed
```
## Extract antifungal SNPs from vcf files
1. Index vcf.gz files with tabix
```
module load bcftools
module load htslib

for file in *.vcf.gz; do
    tabix -f -p vcf "${file}.gz"
done
```
2. create script to merge vcf files
```
vi merge_vcf.sh
```
-Type I to edit and paste:
```
#!/bin/bash
module load bcftools
module load htslib

# Step 1: Extract resistance gene SNPs from each .vcf.gz
for file in *.vcf.gz; do
    base=$(basename "$file" .vcf.gz)
    echo "Extracting SNPs from $file..."
    bcftools view -R resistance_genes.bed "$file" -Oz -o "${base}.resistance.vcf.gz"
    tabix -p vcf "${base}.resistance.vcf.gz"
done

# Step 2: Merge all filtered VCFs
echo "Merging all resistance VCFs..."
bcftools merge *.resistance.vcf.gz -Oz -o merged_resistance.vcf.gz
tabix -p vcf merged_resistance.vcf.gz

# Step 3: Convert merged VCF to SNP matrix
echo "Generating SNP matrix..."
bcftools query -f '%CHROM\t%POS[\t%GT]\n' merged_resistance.vcf.gz > snp_matrix.tsv

echo "All done!"
```
- run script:
```
bash merge_vcf.sh
```
# Let convert vcf to phyllip to run phylogenetic analysis (avoids mafft)
```
module load python/3.8.6 
python vcf2phylip.py -i merged_resistance.vcf.gz -m2
```
- NOTE that the python script vcf2phylip.py has been added to this repository and will be downloaded to the HPC when you clone it.


# Run Raxml
```
module load RAxML
raxmlHPC -s merged_resistance.min2.phy -n antifungal -m GTRCAT -p 12345 -b 12345 -# 100 -f a
```
- change merged_resistance.min2.phy with name of file produced by python script above

# Visualize tree in iTOL
https://itol.embl.de/upload.cgi

- Make sure you display bootstrap support value in nodes.

**Output Files**
- After RAxML runs, you’ll get files like:

- RAxML_bestTree.candida_tree — best-scoring ML tree

- RAxML_bootstrap.candida_tree — bootstrap replicates

- RAxML_bipartitions.candida_tree — best tree with bootstrap support


