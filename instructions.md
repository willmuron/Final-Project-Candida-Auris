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
Will, copy and past script here
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
