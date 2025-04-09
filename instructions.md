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
bcftools mpileup -Ou -f reference.fasta ${SAMPLE}.bam | bcftools call -mv -Ob -o ${SAMPLE}.bcf

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

# Use mafft to align
   - Assume you've collected all target gene sequences from all samples into a    single FASTA file:

-   Copy code

```
antifungal_genes.fasta
```

- Create a script called run_mafft.sh:

```
#!/bin/bash
#SBATCH --job-name=mafft_alignment
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/mafft.out
#SBATCH --mail-user=your_email@domain.edu
#SBATCH --mail-type=END,FAIL

module load mafft/7.475

echo "Starting MAFFT alignment..."

mafft --thread 8 --auto antifungal_genes.fasta > antifungal_genes_aligned.fasta

echo "Alignment completed."
```

# Use raxML to build phylogeny
