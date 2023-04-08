For a repducible run the count matrices can be obtained from ncbi GEO using:

```sh
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE218nnn/GSE225nnn/suppl/GSE225599_RAW.tar
```
Then upacked with:
```sh
tar -xf GSE225599_RAW.tar
```

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


The raw fastq files obtained from SRA using:
```sh
prefetch -v SRR------
fastq-dump --outdir [dst] --split-files [dest]
```
