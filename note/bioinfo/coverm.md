# coverm

[github - wwood/CoverM](https://github.com/wwood/CoverM)
[coverm contig usage](https://wwood.github.io/CoverM/coverm-contig.html)
[coverm genome usage](https://wwood.github.io/CoverM/coverm-genome.html)

Reads per base = Read Count / Length
RPKM(read per kilobase milion) = (Reads per base * 1000) * 1000000 / sum(Count) 
TPM(transcripts per million) マニュアルには`rpkm/total_of_rpkm * 10^6` = (Reads per base * 1000) * 1000000 / sum((Reads per base * 1000) とあるけどなんか合わない。

複数サンプルもいっぺんにかけられる。（contigのLengthは同じなのにサンプル毎に出るけど）
出力のtsvのヘッダーはfastqファイルパスがそのままついてややこしいので、後で変えておくとキレイ。

```
#!/bin/sh
#SBATCH -c 22
#SBATCH -p c6420

source ~/miniconda3/bin/activate bio

mkdir -m 777 assembly/coverm_tmp
TMPDIR=assembly/coverm_tmp

coverm contig \
  --threads 22 \
  --mapper bwa-mem \
  --min-read-percent-identity 0.90 \
  --min-read-aligned-percent  0.50 \
  --exclude-supplementary \
  --methods length count covered_fraction reads_per_base rpkm tpm \
  --coupled \
   fastq/C01_qc_hg38unmap_R1.fastq.gz fastq/C01_qc_hg38unmap_R2.fastq.gz \
   fastq/C03_qc_hg38unmap_R1.fastq.gz fastq/C03_qc_hg38unmap_R2.fastq.gz \
  --reference   assembly/my_catalogue_orf.fasta.gz \
  --output-file assembly/coverm_tmp.tsv

rm -r assembly/coverm_tmp

cat assembly/coverm_tmp.tsv | sed "1s/my_catalogue_orf.fasta.gz\///g" | sed "1s/_qc_hg38unmap_R1.fastq.gz//g" | sed "1s/ /_/g" >assembly/coverm.tsv
rm  assembly/coverm_tmp.tsv
```
