# itsexpress

QIIMEの仲間みたい。

[USDA-ARS-GBRU/itsxpress](https://github.com/usda-ars-gbru/itsxpress)
[Q2-ITSxpress: A tutorial on a QIIME 2 plugin to trim ITS sequences](https://forum.qiime2.org/t/q2-itsxpress-a-tutorial-on-a-qiime-2-plugin-to-trim-its-sequences/5780)


## インストール

```
conda install -c bioconda itsxpress
```

↑で大概は動いたけど、ペアエンドでかけて出力をマージしたリードを出力しようとした場合に、エラーになった。

```
# エラーメッセージ
07-07 09:57 root         ERROR    ITSxpress terminated with errors. See the log file for details.
07-07 09:57 root         ERROR    Command '['bbmerge.sh', 'in=IMET135_cut_R1.fastq.gz', 'in2=IMET135_cut_R2.fastq.gz', 'out=/tmp/itsxpress_5aob4ynb/seq.fq.gz', 't=22', 'maxmismatches=40', 'maxratio=0.3']' returned non-zero exit status 1.
```

bbmapのバージョンによってエラーが出みたいで、たくさん出ているitsxpressとその当時のbbmapを選ぶと動いた↓

```
conda create -n itsxpress -c bioconda itsxpress=1.8.0=py_1
conda install -c bioconda bbmap=38.86=h1296035_0
```

quiime2を入れてそこにitsxpressを足す場合？

```
qiime2-2022.2：https://docs.qiime2.org/2022.2/install/native/#miniconda
wget https://data.qiime2.org/distro/core/qiime2-2022.2-py38-linux-conda.yml
conda env create -n qiime2-2022.2 --file qiime2-2022.2-py38-linux-conda.yml
conda install -c bioconda itsxpress
qiime dev refresh-cache
```

## 動作メモ！

- ペアの場合でR1, R2それぞれで出力する場合、R1（またはR2）がITS全長をまたいでいても後半部分はトリムしてくれない！（ITS1を切り出した場合、R1は18S側のみ、R2は5.8S側のみのトリムとなる）
- 両端が見つかった場合のみITSとみなされる？シングルの場合、またいでいれば両方ともカットしてくれた。またいでない場合（一方の端しか見つからない場合）は抽出してもらえないみたい。。（未確認だけど、ペアでも片方の端しか見つからない場合は抽出してくれない仕組みかも）
- R2だけでシングルでかけるとまたいでいてもITS部分を把握してもらえない？+ストランドだけ探されるみたい。


## ヘルプメッセージ

```
>itsxpress -h
usage: itsxpress [-h] --fastq FASTQ [--single_end] [--fastq2 FASTQ2] --outfile OUTFILE
                 [--outfile2 OUTFILE2] [--tempdir TEMPDIR] [--keeptemp] --region {ITS2,ITS1,ALL}
                 [--taxa {Alveolata,Bryophyta,Bacillariophyta,Amoebozoa,Euglenozoa,Fungi,Chlorophyta,Rhodophyta,Phaeophyceae,Marchantiophyta,Metazoa,Oomycota,Haptophyceae,Raphidophyceae, Rhizaria,Synurophyceae,Tracheophyta,Eustigmatophyceae,All}]
                 [--cluster_id CLUSTER_ID] [--reversed_primers] [--log LOG] [--threads THREADS]

ITSxpress: A python module to rapidly trim ITS amplicon sequences from Fastq files.

optional arguments:
  -h, --help            show this help message and exit
  --fastq FASTQ, -f FASTQ
                        A .fastq, .fq, .fastq.gz or .fq.gz file. Interleaved or not.
  --single_end, -s      A flag to specify that the FASTQ file is single-ended (not paired).
                        Default is false.
  --fastq2 FASTQ2, -f2 FASTQ2
                        A .fastq, .fq, .fastq.gz or .fq.gz file. representing read 2 (optional)
  --outfile OUTFILE, -o OUTFILE
                        the trimmed Fastq file, if it ends in 'gz' it will be gzipped
  --outfile2 OUTFILE2, -o2 OUTFILE2
                        the trimmed read 2 Fastq file, if it ends in 'gz' it will be gzipped. If
                        provided, reads will be returned unmerged.
  --tempdir TEMPDIR     The temp file directory
  --keeptemp            Should intermediate files be kept?
  --region {ITS2,ITS1,ALL}
  --taxa {Alveolata,Bryophyta,Bacillariophyta,Amoebozoa,Euglenozoa,Fungi,Chlorophyta,Rhodophyta,Phaeophyceae,Marchantiophyta,Metazoa,Oomycota,Haptophyceae,Raphidophyceae, Rhizaria,Synurophyceae,Tracheophyta,Eustigmatophyceae,All}
                        The taxonomic group sequenced.
  --cluster_id CLUSTER_ID
                        The percent identity for clustering reads range [0.99-1.0], set to 1 for
                        exact dereplication.
  --reversed_primers, -rp
                        Primers are in reverse orientation as in Taylor et al. 2016,
                        DOI:10.1128/AEM.02576-16. If selected ITSxpress returns trimmed reads
                        flipped to the forward orientation
  --log LOG             Log file
  --threads THREADS     Number of processor threads to use.
```

## example

```
# シングル
itsxpress \
 --region ITS1 \
 --taxa Fungi \
 --threads 22 \
 --single_end \
 --fastq    IMET135_assembled.fastq.gz \
 --fastq2   IMET135_cut_R2.fastq.gz \
 --outfile  IMET135_assembled_its1.fastq.gz \
 --log      IMET135_assembled_its1.log
 
# ペア（R1、R2別々に出力） 
itsxpress \
 --region ITS1 \
 --taxa Fungi \
 --threads 22 \
 --fastq    IMET135_cut_R1.fastq.gz \
 --fastq2   IMET135_cut_R2.fastq.gz \
 --outfile  IMET135_cut_R1_its1.fastq.gz \
 --outfile2 IMET135_cut_R2_its1.fastq.gz \
 --log      IMET135_cut_R1R2_its1.log

# ペア（R1、R2マージで出力） 
itsxpress \
 --region ITS1 \
 --taxa Fungi \
 --threads 22 \
 --fastq    IMET135_cut_R1.fastq.gz \
 --fastq2   IMET135_cut_R2.fastq.gz \
 --outfile  IMET135_cut_R1R2merged_its1.fastq.gz \
 --log      IMET135_cut_R1R2merged_its1.log
```