# MEGAHIT

[GitHub](https://github.com/voutcn/megahit)

```sh
conda install -c bioconda megahit
```

## 使い方

入力fastqは.gz可。アセンブル結果は`outdir_name/final.contigs.fa`として出力される。

```sh
megahit \
  --memory 4 \
  --num-cpu-threads 22 \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -o outdir_name \

# 必要に応じて
gzip outdir_name/final.contigs.fa
```