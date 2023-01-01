# CD-HIT

[bioinformatics/CD-HIT](https://bi.biopapyrus.jp/seq/alignment/software/cd-hit.html)

```
conda install -c bioconda cd-hit
```

## 使い方

アミノ酸配列の時は`cd-hit`、塩基配列の時は`cd-hit-est`。
```
cd-hit-est \
  -T 0 \
  -M 0 \
  -c 0.95 \
  -aS 0.9 \
  -g 1 \
  -d 0 \
  -i sample.fasta \
  -o sample_cdhit.fasta

# 必要に応じて
gzip sample_cdhit.fasta
```

ペアエンドの場合（'-P 1'を忘れずに！）
```
cd-hit-est \
  -T 0 \
  -M 0 \
  -c 0.95 \
  -aS 0.9 \
  -g 1 \
  -d 0 \
  -P 1 \
  -i sample.R1.fasta \
  -j sample.R2.fasta \
  -o sample_cdhit.R1.fasta \
  -op sample_cdhit.R2.fasta 
```

主なオプション
|    |    |
| ----  | ---- |
|  -T   |  使用スレッド数  |
|  -M   |  使用メモリ量  |
|  -P   |  ペアエンドかどうか（0:シングル、1:ペア）  |
|  -c   |  クラスタリングするときの一致度  |
|  -n   | word_length, default 10, see user's guide for choosing it |
|  -aS  |  クラスタリングする時のオーバーラップ（"短い方"のどのくらいオーバーラップしているか） |
|  -g   |  TD  |
|  -d   |  出力ファイルの配列名の長さ  |
|  -i   |  入力ファイル  |
|  -j   |  ペアエンドの場合R2の入力ファイル  |
|  -o   |  出力ファイル  |
|  -op  |  ペアエンドの場合R2の出力ファイル  |
| -mask |  特定の文字をマスク（ex. -mask N） |


clstrファイルからリード数を数えるときは、awkで変形するのが速いよ
```
cat all_R1.fastq.clstr | awk -F "\t" '{if($1 ~ /^>/) {current_name=$1} else {print current_name"\t"$0}}' | head
```

cd-hit-est -maskのテスト。Nが入っていてどちらにでも分類できる場合は、先に出た方？に入れられるみたい。
```
cat test.fa
>1
ATGCATGCATGCATGCATGC
>2
ATGCATGCATGCATGNATGN
>3
ATGCATGCATGCATGGATGG

cd-hit-est -mask N -i test.fa -o test_cdhit.fa

cat test_cdhit.fa.clstr
>Cluster 0
0	20nt, >1... *
1	20nt, >2... at +/100.00%
>Cluster 1
0	20nt, >3... *
```


ヘルプ
```sh
cd-hit-est  -h
		====== CD-HIT version 4.8.1 (built on Apr  7 2021) ======

Usage: cd-hit-est [Options] 

Options

   -i	input filename in fasta format, required, can be in .gz format
   -j	input filename in fasta/fastq format for R2 reads if input are paired end (PE) files
 	 -i R1.fq -j R2.fq -o output_R1 -op output_R2 or
 	 -i R1.fa -j R2.fa -o output_R1 -op output_R2 
   -o	output filename, required
   -op	output filename for R2 reads if input are paired end (PE) files
   -c	sequence identity threshold, default 0.9
 	this is the default cd-hit's "global sequence identity" calculated as:
 	number of identical amino acids or bases in alignment
 	divided by the full length of the shorter sequence
   -G	use global sequence identity, default 1
 	if set to 0, then use local sequence identity, calculated as :
 	number of identical amino acids or bases in alignment
 	divided by the length of the alignment
 	NOTE!!! don't use -G 0 unless you use alignment coverage controls
 	see options -aL, -AL, -aS, -AS
   -b	band_width of alignment, default 20
   -M	memory limit (in MB) for the program, default 800; 0 for unlimitted;
   -T	number of threads, default 1; with 0, all CPUs will be used
   -n	word_length, default 10, see user's guide for choosing it
   -l	length of throw_away_sequences, default 10
   -d	length of description in .clstr file, default 20
 	if set to 0, it takes the fasta defline and stops at first space
   -s	length difference cutoff, default 0.0
 	if set to 0.9, the shorter sequences need to be
 	at least 90% length of the representative of the cluster
   -S	length difference cutoff in amino acid, default 999999
 	if set to 60, the length difference between the shorter sequences
 	and the representative of the cluster can not be bigger than 60
   -aL	alignment coverage for the longer sequence, default 0.0
 	if set to 0.9, the alignment must covers 90% of the sequence
   -AL	alignment coverage control for the longer sequence, default 99999999
 	if set to 60, and the length of the sequence is 400,
 	then the alignment must be >= 340 (400-60) residues
   -aS	alignment coverage for the shorter sequence, default 0.0
 	if set to 0.9, the alignment must covers 90% of the sequence
   -AS	alignment coverage control for the shorter sequence, default 99999999
 	if set to 60, and the length of the sequence is 400,
 	then the alignment must be >= 340 (400-60) residues
   -A	minimal alignment coverage control for the both sequences, default 0
 	alignment must cover >= this value for both sequences 
   -uL	maximum unmatched percentage for the longer sequence, default 1.0
 	if set to 0.1, the unmatched region (excluding leading and tailing gaps)
 	must not be more than 10% of the sequence
   -uS	maximum unmatched percentage for the shorter sequence, default 1.0
 	if set to 0.1, the unmatched region (excluding leading and tailing gaps)
 	must not be more than 10% of the sequence
   -U	maximum unmatched length, default 99999999
 	if set to 10, the unmatched region (excluding leading and tailing gaps)
 	must not be more than 10 bases
   -B	1 or 0, default 0, by default, sequences are stored in RAM
 	if set to 1, sequence are stored on hard drive
 	!! No longer supported !!
   -P	input paired end (PE) reads, default 0, single file
 	if set to 1, please use -i R1 -j R2 to input both PE files
   -cx	length to keep after trimming the tail of sequence, default 0, not trimming
 	if set to 50, the program only uses the first 50 letters of input sequence
   -cy	length to keep after trimming the tail of R2 sequence, default 0, not trimming
 	if set to 50, the program only uses the first 50 letters of input R2 sequence
 	e.g. -cx 100 -cy 80 for paired end reads
   -ap	alignment position constrains,  default 0, no constrain
 	if set to 1, the program will force sequences to align at beginings
 	when set to 1, the program only does +/+ alignment
   -p	1 or 0, default 0
 	if set to 1, print alignment overlap in .clstr file
   -g	1 or 0, default 0
 	by cd-hit's default algorithm, a sequence is clustered to the first 
 	cluster that meet the threshold (fast cluster). If set to 1, the program
 	will cluster it into the most similar cluster that meet the threshold
 	(accurate but slow mode)
 	but either 1 or 0 won't change the representatives of final clusters
   -r	1 or 0, default 1, by default do both +/+ & +/- alignments
 	if set to 0, only +/+ strand alignment
   -mask	masking letters (e.g. -mask NX, to mask out both 'N' and 'X')
   -match	matching score, default 2 (1 for T-U and N-N)
   -mismatch	mismatching score, default -2
   -gap	gap opening score, default -6
   -gap-ext	gap extension score, default -1
   -bak	write backup cluster file (1 or 0, default 0)
   -sc	sort clusters by size (number of sequences), default 0, output clusters by decreasing length
 	if set to 1, output clusters by decreasing size
   -sf	sort fasta/fastq by cluster size (number of sequences), default 0, no sorting
 	if set to 1, output sequences by decreasing cluster size
 	this can be very slow if the input is in .gz format
   -h	print this help

   Questions, bugs, contact Weizhong Li at liwz@sdsc.edu
   For updated versions and information, please visit: http://cd-hit.org
                                                    or https://github.com/weizhongli/cdhit

   cd-hit web server is also available from http://cd-hit.org

   If you find cd-hit useful, please kindly cite:

   "CD-HIT: a fast program for clustering and comparing large sets of protein or nucleotide sequences", Weizhong Li & Adam Godzik. Bioinformatics, (2006) 22:1658-1659
   "CD-HIT: accelerated for clustering the next generation sequencing data", Limin Fu, Beifang Niu, Zhengwei Zhu, Sitao Wu & Weizhong Li. Bioinformatics, (2012) 28:3150-3152

```
