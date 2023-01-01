# SeqKit

[SeqKit - Usage](https://bioinf.shenwei.me/seqkit/usage/)
[みんなのためのバイオインフォマティクス - 
fasta, fastq ファイルを扱うときに便利なseqkit](https://bioinfo-nanihitotsu.weebly.com/software/category/seqkit)

```
conda install -c bioconda megahit
```

ツール一覧
```
>seqkit help
SeqKit -- a cross-platform and ultrafast toolkit for FASTA/Q file manipulation

Version: 0.16.1

Author: Wei Shen <shenwei356@gmail.com>

Documents  : http://bioinf.shenwei.me/seqkit
Source code: https://github.com/shenwei356/seqkit
Please cite: https://doi.org/10.1371/journal.pone.0163962

Usage:
  seqkit [command]

Available Commands:
  amplicon        retrieve amplicon (or specific region around it) via primer(s)
  bam             monitoring and online histograms of BAM record features
  common          find common sequences of multiple files by id/name/sequence
  concat          concatenate sequences with same ID from multiple files
  convert         convert FASTQ quality encoding between Sanger, Solexa and Illumina
  duplicate       duplicate sequences N times
  faidx           create FASTA index file and extract subsequence
  fish            look for short sequences in larger sequences using local alignment
  fq2fa           convert FASTQ to FASTA
  fx2tab          convert FASTA/Q to tabular format (with length/GC content/GC skew)
  genautocomplete generate shell autocompletion script (bash|zsh|fish|powershell)
  grep            search sequences by ID/name/sequence/sequence motifs, mismatch allowed
  head            print first N FASTA/Q records
  head-genome     print sequences of the first genome with common prefixes in name
  help            Help about any command
  locate          locate subsequences/motifs, mismatch allowed
  mutate          edit sequence (point mutation, insertion, deletion)
  pair            match up paired-end reads from two fastq files
  range           print FASTA/Q records in a range (start:end)
  rename          rename duplicated IDs
  replace         replace name/sequence by regular expression
  restart         reset start position for circular genome
  rmdup           remove duplicated sequences by id/name/sequence
  sample          sample sequences by number or proportion
  sana            sanitize broken single line fastq files
  scat            real time recursive concatenation and streaming of fastx files
  seq             transform sequences (revserse, complement, extract ID...)
  shuffle         shuffle sequences
  sliding         sliding sequences, circular genome supported
  sort            sort sequences by id/name/sequence/length
  split           split sequences into files by id/seq region/size/parts (mainly for FASTA)
  split2          split sequences into files by size/parts (FASTA, PE/SE FASTQ)
  stats           simple statistics of FASTA/Q files
  subseq          get subsequences by region/gtf/bed, including flanking sequences
  tab2fx          convert tabular format to FASTA/Q format
  translate       translate DNA/RNA to protein sequence (supporting ambiguous bases)
  version         print version information and check for update
  watch           monitoring and online histograms of sequence features

Flags:
      --alphabet-guess-seq-length int   length of sequence prefix of the first FASTA record based on which seqkit guesses the sequence type (0 for whole seq) (default 10000)
  -h, --help                            help for seqkit
      --id-ncbi                         FASTA head is NCBI-style, e.g. >gi|110645304|ref|NC_002516.2| Pseud...
      --id-regexp string                regular expression for parsing ID (default "^(\\S+)\\s?")
      --infile-list string              file of input files list (one file per line), if given, they are appended to files from cli arguments
  -w, --line-width int                  line width when outputing FASTA format (0 for no wrap) (default 60)
  -o, --out-file string                 out file ("-" for stdout, suffix .gz for gzipped out) (default "-")
      --quiet                           be quiet and do not show extra information
  -t, --seq-type string                 sequence type (dna|rna|protein|unlimit|auto) (for auto, it automatically detect by the first sequence) (default "auto")
  -j, --threads int                     number of CPUs. (default value: 1 for single-CPU PC, 2 for others. can also set with environment variable SEQKIT_THREADS) (default 2)

Use "seqkit [command] --help" for more information about a command.
```

## 1番使うseqkit stats

`-T`でタブ区切り

```

```

## 配列の抽出(grep)

配列名から抽出

```
seqkit grep -nrp "M00242:626:000000000-G6P2D:1:1104:16846:28159" sample.fastq.gz
```

一致「しないもの」を抽出するときは-vオプション

```
seqkit grep -v -nrp "KT715474.1" ntF-ITS1.fasta > ntF-ITS1_modified.fasta
```

配列名のリストを使って抽出（便利！）

```
head list.txt
read1
read5
reaf6

seqkit grep -n -f list.txt sample.fasta.gz > sample_grep.fasta
```

配列から検索

```
# ギャップを含むリードを抽出
seqkit grep -s -p - sample.fasta.gz

# アダプター配列があるか確認
seqkit grep -s -p CTGTCTCTTATACACATCTCCGAGCCCACGAGAC NKTAall_R1.fastq.gz | seqkit stats
file  format  type  num_seqs     sum_len  min_len  avg_len  max_len
-     FASTQ   DNA    153,714  38,582,214      251      251      251
```

## ランダム抽出

```
# 割合で
seqkit sample -p 0.1 seqs.fq.gz

# リード数で
seqkit sample -n 1000 seqs.fq.gz

# ファイルが大きいとき、-nを使わない（全部読み込むため時間がかかる）。代わりに-pのあとseqkit headを使う
seqkit sample -p 0.1 seqs.fq.gz | seqkit head -n N
```


## アミノ酸配列に変換(translate)
```
seqkit translate --frame 1 --transl-table 11 input.fasta.gz > input_aa.fa
```

## リードをランダムサンプリング(samlple)
```
seqkit sample -p 0.1 input.fastq > output.fastq

# 主なオプション
-n, --number int         sample by number (result may not exactly match), DO NOT use on large FASTQ files.
-p, --proportion float   sample by proportion
-s, --rand-seed int      rand seed (default 11) # ペアリードはrandom seedを同じ値にする。
-j, --threads int        number of CPUs. (default value: 1 for single-CPU PC, 2 for others. can also set with environment variable SEQKIT_THREADS) (default 2)
```

## reverse complementary (相補鎖にしてさらにリバースにする）

```
seqkit seq --seq-type dna -pr input.fa > output.fa
```

## 配列名や配列の置換（seqkit replace）

```
# 配列名冒頭にサンプル名を追加
for i in $(cat ../sample.txt)
cat ${i}.R1.fastq \
| seqkit replace -p '^' -r ${i}_ >${i}_rename.R1.fastq
```

```
# 配列名にサンプル名を追加して１ファイルにまとめるとき
for i in $(cat sample.txt)
do zcat 1_fastq/${i}_R1.fastq.gz | seqkit replace -p '^' -r ${i}_ >>sppFoundIn07/NKTAall_R1.fastq
   zcat 1_fastq/${i}_R2.fastq.gz | seqkit replace -p '^' -r ${i}_ >>sppFoundIn07/NKTAall_R2.fastq
done
gzip sppFoundIn07/NKTAall_R1.fastq
gzip sppFoundIn07/NKTAall_R2.fastq
```

## 配列/配列名の重複を除去（seqkit rmdup）

```
# 配列名の重複を除去（実行字に重複数を教えてくれるので、重複がないことを確かめるときもよい）
seqkit rmdup -n all_Candida_albicans_NR_125332_mapped.R1.fastq_1 -o tmp.fastq
```

## seqkit amplicon
[Trim fastq after and before motif occurance](https://www.biostars.org/p/397146/)

特定のモチーフ（プライマー）に挟まれた領域内を取り出したりする。モチーフの前後何ベースの取り出しなどもできる。一致する箇所がないリードは捨てられる。ミスマッチ数は設定できる。-Rの方はreverce complementで書くこと。
```
cat test2.fa
>1
AAATATGCATGCGCCC
>2
AAAAATATGCATGCGCCCCC

seqkit amplicon -F AAAT -R GGGC test2.fa
[INFO] 1 primer pair loaded
>1
AAATATGCATGCGCCC
>2
AAATATGCATGCGCCC

seqkit amplicon -r 5:-5 -F AAAT -R GGGC test2.fa
[INFO] 1 primer pair loaded
>1
ATGCATGC
>2
ATGCATGC


seqkit amplicon -F AAAAAT -R GGGGGC test2.fa
[INFO] 1 primer pair loaded
>2
AAAAATATGCATGCGCCCCC

# デフォではreverce complementも探される。出力はreverce complementに変換されている。
seqkit amplicon -F GGGC -R AAAT test2.fa
[INFO] 1 primer pair loaded
>1
GGGCGCATGCATATTT
>2
GGGCGCATGCATATTT

# reverce complementは探したくない時
seqkit amplicon -F GGGC -R AAAT --only-positive-strand test2.fa
[INFO] 1 primer pair loaded

# 片方だけだとその配列だけ出る
seqkit amplicon -F AAAT  test2.fa
[INFO] 1 primer pair loaded
>1
AAAT
>2
AAAT

# 片方だけ指定してそこを含めた後ろ10baseを取り出す場合
seqkit amplicon -F AAAT -r 1:10 test2.fa 
[INFO] 1 primer pair loaded
>1
AAATATGCAT
>2
AAATATGCAT
```

両端（プライマー）配列はタブ区切りファイルでも指定できる（複数可）。
```
cat test2_primer.tsv
p1	AAAT	GGGC

seqkit amplicon --primer-file test2_primer.tsv test2.fa
[INFO] 1 primer pair loaded
>1
AAATATGCATGCGCCC
>2
AAATATGCATGCGCCC
```

bedファイルにする
```
seqkit amplicon -F AAAT -R GGGC --bed  test2.fa
[INFO] 1 primer pair loaded
1	0	16	.	0	+	AAATATGCATGCGCCC
2	2	18	.	0	+	AAATATGCATGCGCCC

seqkit amplicon -F AAAT  --bed  test2.fa
[INFO] 1 primer pair loaded
1	0	4	.	0	+	AAAT
2	2	6	.	0	+	AAAT
```