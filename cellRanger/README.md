Welcome. The instructions provided here are designed to help you create a reference then run raw (fastq) single cell sequencing data through a cellranger pipeline.

# Steps to create reference genome and run cellranger:
1. [Install cellranger](#1-install-cellranger)
2. [Download and prepare a reference genome](#2-download-and-prepare-a-reference-genome)
3. [Download and prepare the GTF annotation file](#3-download-and-prepare-the-GTF-files)
4. [Convert the GTF and genome to a cellranger reference](#4-convert-gtf-file-and-genome-to-cellranger-reference-file)
5. [Run cellranger counts to align data](#6-run-cellranger-counts)

## 1. Install cellranger

Unpack the tar ball after it finishes downloading.
``` sh
tar -xzvf cellranger*tar.gz
```
#### Test that you have access to cellranger:
```sh
export PATH=/projects/$USER/software/cellranger-7.0.0:$PATH
cellranger
```

## 2. Download and prepare a reference genome:

#### Get the reference files:
```sh
rsync -avzP rsync://ftp.ensembl.org/ensembl/pub/release-104/fasta/canis_lupus_familiaris/dna/*.dna.toplevel*.fa.gz .
gunzip *.dna.toplevel*.fa.gz
```		
## 3. Download and prepare the GTF files:

#### Pull the GTF from ensembl (current as of Dec. 16, 2021):
```sh
rsync -avzP rsync://ftp.ensembl.org/ensembl/pub/release-104/gtf/canis_lupus_familiaris/*.gtf.gz  .
gunzip *.gtf.gz
```
#### Filter the GTF file with cellranger mkgtf:

```sh
bash mkgtf.sh > mkgtf.log 2>&1 &
```

## 4. Convert the gtf file and genome to cellranger reference file:

```sh
sbatch cute_cellrngr_mkref.sbatch
```	

## 5. Run cellranger counts

```sh
sbatch cute_cellrngrCnts.sbatch
```

