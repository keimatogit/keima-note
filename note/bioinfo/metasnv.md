# metaSNV

[github](https://github.com/metasnv-tool/metaSNV)
[metaSNV](https://metasnv.embl.de/)
[manual](https://github.com/metasnv-tool/metaSNV/blob/master/documentation/metaSNV_v2_manual.pdf)


## インストールなど

```
conda create -n metasnv -c bioconda metasnv
```

ちゃんとインストールされてることを確認

```
metaSNV.py --help
metaSNV_Filtering.py --help
metaSNV_DistDiv.py --help
metaSNV_subpopr.R --help
```

レファレンスはユーザーで設定できるが、１種につき１ゲノムのものを使用すること。
ProGenomes2レファレンスとが用意されているのでダウンロード。（1種につき1つの代表ゲノムのサブセット（代表ゲノムのうち最も長いもの）で、２つで25 GB程度とのこと。）

```
wget http://swifter.embl.de/~ralves/metaSNV_reference_data/progenomes2_speciesReps_genomes.fna
wget http://swifter.embl.de/~ralves/metaSNV_reference_data/progenomes2_speciesReps_annotations.txt
```

bwaのインデックス作成（マッピングソフトはなんでも良いが、ここではbwaを使用）

```
bwa index -p bwa_index/progenomes2_speciesReps progenomes2_speciesReps_genomes.fna
```


## ヘルプ

metaSNV.py --help

```
usage: metaSNV.py [-h] [--db_ann DB_ANN_FILE] [--print-commands] [--threads INT] [--n_splits INT]
                  [--use_prev_cov] [--min_pos_cov INT] [--min_pos_snvs INT]
                  DIR FILE REF_DB_FILE

Compute SNV profiles

positional arguments:
  DIR                   The output directory that metaSNV will create e.g. "outputs". Can be a
                        path.
  FILE                  File with an input list of bam files, one file per line
  REF_DB_FILE           reference multi-sequence FASTA file used for the alignments.

options:
  -h, --help            show this help message and exit
  --db_ann DB_ANN_FILE  Database gene annotation.
  --print-commands      Instead of executing the commands, simply print them out
  --threads INT         Number of jobs to run simmultaneously. Will create same number of splits,
                        unless n_splits set differently.
  --n_splits INT        Number of bins to split ref into
  --use_prev_cov        Use "cov/" and "outputs.all_cov.tab" and "outputs.all_perc.tab" data
                        produced by previous metaSNV run
  --min_pos_cov INT     minimum coverage (mapped reads) per position for snpCall.
  --min_pos_snvs INT    minimum number of non-reference nucleotides per position for snpCall.
```

## Workflow


### マッピング

レファレンスにマッピングしてbamファイルを作成する（bamファイルはソートしておくこと）。マッピングツールはなんでも良い（BWA, bowtie, Minimap2など）が、ngslessでBWAを使う方法がマニュアルで紹介されていた。ここではbwaを使用。 

BWA-mem

```

```

BAMファイルができたら、ファイル名のリストを作成する。

bam_list.txt

```
/path/to/sample1.sort.bam
/path/to/sample2.sort.bam
/path/to/sample3.sort.bam
```

Call SNVs

```
metaSNV.py output_dir/ bam_list.txt /path/to/reference.fna [options]
```

SNV Post-Processing: Filtering & Analysis

```
metaSNV_Filtering.py output_dir [options]
metaSNV_DistDiv.py --filt output_dir/filtered/pop [options]
```

Subpopulation detection

```
metaSNV_subpopr.R -i output_dir [options]
```