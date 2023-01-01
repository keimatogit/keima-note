# Trinity

- [Trinity Wiki Home](https://github.com/trinityrnaseq/trinityrnaseq/wiki)
- [Trinity Wiki - Trinity Transcript Quantification](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Transcript-Quantification)

インストール
ちょっと手間取ったみたい。専用環境を作った方がいいと思う。
```
conda acreate -n trinity
conda activate trinity
conda install -c bioconda bowtie samtools trinity=2.8.5=h8b12597_5 # 全部いっぺんにインストールしたら動いた。
conda install -c bioconda rsem bowtie2
```


実行shの例
```
#!/bin/bash
#SBATCH -c 88

source ~/miniconda3/bin/activate trinity

# seqType: [fq, fa]
# output must include "trinity"
# more than one fasta/fastq file: separated by commas, no spaces (R1=1.fastq\,2_fastq)
# assebled reads are in the outpdir/Trinity.fasta
# normaly, --max_memory 50G, --CPU 22 (-c 22) ??
# making temporary file (*.readcount) in the same place of R1 and R2 (will be removed)

Trinity \
 --seqType fq \
 --max_memory 150G \
 --CPU    22 \
 --left \
0_fastq/CD11_clean_R1.fastq.gz,\
0_fastq/CD13_clean_R1.fastq.gz,\
0_fastq/CD14_clean_R1.fastq.gz,\
0_fastq/CD15_clean_R1.fastq.gz,\
0_fastq/CD16_clean_R1.fastq.gz,\
0_fastq/CD17_clean_R1.fastq.gz,\
0_fastq/CD1_clean_R1.fastq.gz,\
0_fastq/CD2_clean_R1.fastq.gz,\
0_fastq/CD3_clean_R1.fastq.gz,\
0_fastq/CD4_clean_R1.fastq.gz,\
0_fastq/CD9_clean_R1.fastq.gz,\
0_fastq/HC10_clean_R1.fastq.gz,\
0_fastq/HC11_clean_R1.fastq.gz,\
0_fastq/HC1_clean_R1.fastq.gz,\
0_fastq/HC2_clean_R1.fastq.gz,\
0_fastq/HC3_clean_R1.fastq.gz,\
0_fastq/HC5_clean_R1.fastq.gz,\
0_fastq/HC7_clean_R1.fastq.gz,\
0_fastq/HC8_clean_R1.fastq.gz \
 --right \
0_fastq/CD11_clean_R2.fastq.gz,\
0_fastq/CD13_clean_R2.fastq.gz,\
0_fastq/CD14_clean_R2.fastq.gz,\
0_fastq/CD15_clean_R2.fastq.gz,\
0_fastq/CD16_clean_R2.fastq.gz,\
0_fastq/CD17_clean_R2.fastq.gz,\
0_fastq/CD1_clean_R2.fastq.gz,\
0_fastq/CD2_clean_R2.fastq.gz,\
0_fastq/CD3_clean_R2.fastq.gz,\
0_fastq/CD4_clean_R2.fastq.gz,\
0_fastq/CD9_clean_R2.fastq.gz,\
0_fastq/HC10_clean_R2.fastq.gz,\
0_fastq/HC11_clean_R2.fastq.gz,\
0_fastq/HC1_clean_R2.fastq.gz,\
0_fastq/HC2_clean_R2.fastq.gz,\
0_fastq/HC3_clean_R2.fastq.gz,\
0_fastq/HC5_clean_R2.fastq.gz,\
0_fastq/HC7_clean_R2.fastq.gz,\
0_fastq/HC8_clean_R2.fastq.gz \
 --output 1_trinity \
 --full_cleanup

```

## Transcript Quantification

220311 F3572でalign_and_estimate_abundance.plとbowtie2->RSEM両方やってみたらちょっとだけ数値が違ったんだけど何。
align_and_estimate_abundance.pl標準出力によると、bowtie2実行時に-X 800オプションも入ってる？からかも。-Xは有効なペアエンドアラインメントの最大フラグメント長でデフォは500。

### align_and_estimate_abundance.pl

- miniconda3/envs/trinity/bin/にある。
- マッピングから通しでできる。マッピングは別にやってbamファイルを渡すこともできる。
- --prep_referenceで作られるbowtie2やrsemのインデックスは--transcriptsと同じ場所できる。
- bowtie2に渡すオプションのデフォルトは、`-no-mixed --no-discordant --gbar 1000 --end-to-end -k 200 -X 800`。
- --samples_fileで4カラムのサンプルファイル一覧（tab-delimited text file indicating biological replicate relationships.）を渡すこともできる（align_and_estimate_abundance.pl -h 参照）


```sh
# レファレンス準備
# （--transcriptsと同じフォルダにTrinity.fasta.bowtie2.*.bt2とTrinity.fasta.RSEM.*が作られる）
# （bowtie2-buildやrsem-prepare-referenceで作られるやつ）
align_and_estimate_abundance.pl \
  --thread_count 4 \
  --transcripts ../1_trinity/Trinity.fasta \
  --est_method RSEM \
  --aln_method bowtie2 \
  --trinity_mode \
  --prep_reference

# 実行
align_and_estimate_abundance.pl \
  --thread_count 4 \
  --transcripts ../1_trinity/Trinity.fasta \
  --seqType fq \
  --left  ../1_sortmerna/lard_1_clean_R1.fastq.gz \
  --right ../1_sortmerna/lard_1_clean_R2.fastq.gz \
  --output_dir 4_align_and_estimate_abundance \
  --est_method RSEM \
  --aln_method bowtie2 \
  --gene_trans_map ../1_trinity/Trinity.fasta.gene_trans_map
```
output_dirに以下のファイルができる。
- bowtie2.bam
- bowtie2.bam.ok
- RSEM.isoforms.results
- RSEM.stat
- bowtie2.bam.for_rsem.bam
- RSEM.genes.results
- RSEM.isoforms.results.ok

サンプルファイルを指定下場合、カレントディレクトリ（指定すればoutput_dir?）にサンプル名のフォルダができて、その中に↑mpファイルができる。

```
#  Example usage
#
#   ## Just prepare the reference for alignment and abundance estimation
#
#   align_and_estimate_abundance.pl --transcripts Trinity.fasta --est_method RSEM --aln_method bowtie --trinity_mode --prep_reference
#
#   ## Run the alignment and abundance estimation (assumes reference has already been prepped, errors-out if prepped reference not located.)
#
#    align_and_estimate_abundance.pl --transcripts Trinity.fasta --seqType fq --left reads_1.fq --right reads_2.fq --est_method RSEM --aln_method bowtie --trinity_mode --output_dir rsem_outdir
#
##  ## prep the reference and run the alignment/estimation
#
#    align_and_estimate_abundance.pl --transcripts Trinity.fasta --seqType fq --left reads_1.fq --right reads_2.fq --est_method RSEM --aln_method bowtie --trinity_mode --prep_reference --output_dir rsem_outdir
#
#   ## Use a samples.txt file:
#
#    align_and_estimate_abundance.pl --transcripts Trinity.fasta --est_method RSEM --aln_method bowtie2 --prep_reference --trinity_mode --samples_file samples.txt --seqType fq  
#

```

`--samples_file`の例
```
CD	CD11	../0_fastq/CD11_clean_R1.fastq.gz	../0_fastq/CD11_clean_R2.fastq.gz
CD	CD13	../0_fastq/CD13_clean_R1.fastq.gz	../0_fastq/CD13_clean_R2.fastq.gz
CD	CD14	../0_fastq/CD14_clean_R1.fastq.gz	../0_fastq/CD14_clean_R2.fastq.gz
HC	HC10	../0_fastq/HC10_clean_R1.fastq.gz	../0_fastq/HC10_clean_R2.fastq.gz
HC	HC11	../0_fastq/HC11_clean_R1.fastq.gz	../0_fastq/HC11_clean_R2.fastq.gz
HC	HC1	../0_fastq/HC1_clean_R1.fastq.gz	../0_fastq/HC1_clean_R2.fastq.gz
```

### bowtie2 ~ RSEM