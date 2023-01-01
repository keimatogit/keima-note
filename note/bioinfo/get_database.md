# get database

## hg38(uscs)

```
cd ~/db/ucsc/hg38

# fasta
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

# genes (refseq)
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz
```


## ヒトトランスクリプトーム（GRCh38, RefSeq Transcripts, NCBI）

[Human Genome Resources at NCBI](https://www.ncbi.nlm.nih.gov/genome/guide/human/)
221214ダウンロード

```
cd db/ncbi/fasta
wget  https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_rna.fna.gz
```

## CAZy

[dbCAN](https://bcb.unl.edu/dbCAN2/)

```
cd ~/db/dbCAN
wget https://bcb.unl.edu/dbCAN2/download/CAZyDB.08062022.fa
wget https://bcb.unl.edu/dbCAN2/download/dbCAN-HMMdb-V11.txt
```