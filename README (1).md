# 🧬 Final Project: Candida auris Data Processing Pipeline

This pipeline automates the process of downloading raw sequencing data from the NCBI Sequence Read Archive (SRA), converting the data to FASTQ format, and running quality control using FastQC.

---

## 📦 Prerequisites

Ensure that the following software and tools are installed or available as modules in your HPC environment:

- [SRA Toolkit v2.10.9](https://github.com/ncbi/sra-tools)
- [FastQC v0.11.9](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

Load them with:

```bash
module load sra-toolkit/2.10.9
module load FastQC/0.11.9
```

---

## 🛠️ Pipeline Overview

This pipeline performs the following steps:

1. **Load Modules** – Ensure tools are available in your environment.
2. **Create Directory Structure** – Set up folders for SRA, FASTQ, and QC files.
3. **Download SRA Files** – Using accession numbers in `accession_list.txt`.
4. **Convert to FASTQ** – Using `fastq-dump` to split and compress reads.
5. **Run FastQC** – Perform quality checks on the FASTQ files.

---

## 📋 Step-by-Step Instructions

### 1. Create an Accession List

Create a file named `accession_list.txt` containing your SRA accession numbers, one per line:

```text
SRR7140069
SRR7140058
SRR7976608
...
```

### 2. Create the Pipeline Script

Create a file named `run_pipeline.sh` and add the following code:

```bash
#!/bin/bash

# Load required modules
module load sra-toolkit/2.10.9
module load FastQC/0.11.9

# Set SRA download directory
export SRA_DIR="/mnt/scratch/${USER}_candida_sra"
mkdir -p "$SRA_DIR"

# Confirm accession list is present
if [ ! -f accession_list.txt ]; then
    echo "❌ accession_list.txt not found!"
    exit 1
fi

# Step 1: Download SRA files
echo "[1] Downloading SRA data..."
prefetch --output-directory "$SRA_DIR" --option-file accession_list.txt

# Step 2: Convert SRA to FASTQ
echo "[2] Converting .sra to .fastq.gz..."
mkdir -p fastq_files
for sra in $(cat accession_list.txt); do
    if [ -f "$SRA_DIR/$sra/$sra.sra" ]; then
        fastq-dump --split-files --gzip --outdir fastq_files "$SRA_DIR/$sra/$sra.sra"
    else
        echo "⚠️ $sra.sra not found — skipping"
    fi
done

# Step 3: Run FastQC
echo "[3] Running FastQC..."
mkdir -p fastqc_reports
fastqc fastq_files/*.fastq.gz --outdir fastqc_reports

echo "✅ Pipeline complete: SRA downloaded, FASTQ files created, QC complete."
```

### 3. Make the Script Executable

```bash
chmod +x run_pipeline.sh
```

### 4. Run the Pipeline

```bash
./run_pipeline.sh
```

---

## 📁 Output

- `fastq_files/`: Contains your `.fastq.gz` files.
- `fastqc_reports/`: Contains HTML and zip files with FastQC quality control reports.
- `/mnt/scratch/$USER_candida_sra/`: Temporary storage of `.sra` files (can be cleaned up later).

---

## 📝 Notes

- Ensure `accession_list.txt` is in the same directory as your script.
- Adjust paths and module versions if needed for your environment.
- Consider deleting `.sra` files after conversion to save storage.

---

## 📚 References

- **SRA Toolkit:** https://github.com/ncbi/sra-tools  
- **FastQC:** https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
