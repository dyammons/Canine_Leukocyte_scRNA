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


To rearrange the file structure for easy loading into Seurat following the code chunks below (Alternatively, manually rearrange):

```sh
touch rearrange.sh
```

Then copy the conents of the below script into the rearrange.sh file
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
Navigate to directory of interest and run for each file you wish to pull down

NOTE: all raw data file are just under 2TB of data
```sh
prefetch -v -max-size=35000000 SRR------
fastq-dump --outdir [dest] --split-files [dest]
```
