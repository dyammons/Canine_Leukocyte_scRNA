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


To reaggrange the file structure for easy loading into R the following bash script is provided:
```sh
while read line
do
    old=$(echo "$line" | cut -f1)
    new=$(echo "$line" | cut -f2)

    barcodes="./$old/barcodes.tsv.gz"
    feats="./$old/features.tsv.gz"
    mtx="./$old/matrix.mtx.gz"

    cp $barcodes ./forGEO_2/$new\_barcodes.tsv.gz
    cp $feats ./forGEO_2/$new\_features.tsv.gz
    cp $mtx ./forGEO_2/$new\_matrix.mtx.gz

done < deCoder.tsv
```
(Alternatively, manually rearrange)


From there you can follow the runScript_PBMC_analysis.R script


### Instructions to obtain raw data from SRA
Navigate to directory of interest and run for each file you wish to pull down

NOTE: all raw data file are just under 2TB of data
```sh
prefetch -v -max-size=35000000 SRR------
fastq-dump --outdir [dest] --split-files [dest]
```
