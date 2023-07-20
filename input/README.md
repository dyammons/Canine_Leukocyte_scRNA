### Instructions to obtain processed data

From terminal (recommended):
```sh
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225599/suppl/GSE225599_bcell.rds.gz #b cells
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225599/suppl/GSE225599_cytotoxic.rds.gz #cd8 t, nk, dn t
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225599/suppl/GSE225599_final_dataSet_H.rds.gz #just healthy
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225599/suppl/GSE225599_final_dataSet_HvO.rds.gz #combined
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225599/suppl/GSE225599_helper.rds.gz #cd4 t
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225599/suppl/GSE225599_myeloid.rds.gz #myeloid cells

gunzip *.rds.gz
```
From R (NOTE: be sure to change `dest` to the path you want to download to!):
```r
utils::download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225599/suppl/GSE225599_bcell.rds.gz", dest = "/pwd/to/dir/bcell.rds.gz")
utils::download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225599/suppl/GSE225599_cytotoxic.rds.gz", dest = "/pwd/to/dir/cytotoxic.rds.gz")
utils::download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225599/suppl/GSE225599_final_dataSet_H.rds.gz", dest = "/pwd/to/dir/final_dataSet_H.rds.gz")
utils::download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225599/suppl/GSE225599_final_dataSet_HvO.rds.gz", dest = "/pwd/to/dir/final_dataSet_HvO.rds.gz")
utils::download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225599/suppl/GSE225599_helper.rds.gz", dest = "/pwd/to/dir/helper.rds.gz")
utils::download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225599/suppl/GSE225599_myeloid.rds.gz", dest = "/pwd/to/dir/myeloid.rds.gz")
```
You can then use the R function `GEOquery::gunzip`, terminal `gunzip`, or prefered method to unzip the files

### Instructions to obtain count matrices from NCBI
Navigate to directory in which you wish to download the data and run the following:
```sh
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225599/suppl/GSE225599_RAW.tar
```
(Alternatively, you can download the zip folder by visiting the [GSE225599](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225599) page)


Then upacked with:
```sh
tar -xf GSE225599_RAW.tar
```
(Or, right click and extract all)


To rearrange the file structure for easy loading into Seurat, follow the code chunks below (Alternatively, manually rearrange):

First, get the provided deCoder.tsv file in the directory you wish to house the input files then:
```sh
touch rearrange.sh
```

Then, copy the contents of the below script into the rearrange.sh file
```sh
#!/usr/bin/env bash

while read line
do
    old=$(echo "$line" | cut -f1)
    new=$(echo "$line" | cut -f2)
    
    mkdir $new
    
    barcodes="./${old}_barcodes.tsv.gz"
    feats="./${old}_features.tsv.gz"
    mtx="./${old}_matrix.mtx.gz"

    mv $barcodes ./$new/$new\_barcodes.tsv.gz
    mv $feats ./$new/$new\_features.tsv.gz
    mv $mtx ./$new/$new\_matrix.mtx.gz

done < deCoder.tsv
```

After transfering the code, run:
```sh
bash rearrange.sh
```

From there you can follow the runScript_PBMC_analysis.R script

### Instructions to obtain raw data from SRA
Navigate to directory of interest and run for each file you wish to pull down. Then with SRA toolkit installed run:

NOTE: all raw data files are just under 2TB of data
```sh
prefetch -v --max-size=55000000 SRR23525324 #smallest file for ex
fastq-dump --gzip --split-files SRR23525324
```
From there you will have to modify the file names to be compatible with the cellranger input (if using cellranger). It expects `[Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz` -- feel free to reach out if you have trouble
