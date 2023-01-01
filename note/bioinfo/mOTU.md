# mOTUs

- [mOTUs](https://motu-tool.org/index.html)
- [mOTUs wiki](https://github.com/motu-tool/mOTUs/wiki)

インストール
```
conda acreate -n trinity
conda activate trinity
conda install -c bioconda bowtie samtools trinity=2.8.5=h8b12597_5 # 全部いっぺんにインストールしたら動いた。
conda install -c bioconda rsem bowtie2
```

## Taxonomic profile

### 1. Profile one sample

> We suggest to use 4 or 8 threads, note that the time needed for profiling is not scaling linearly with the number of threads.
> Note that the name that is specified with -n will be used in the merge function of motus.
> type of read counts (-y)... The values with insert.* counts the number of inserts (reads) that map to the gene, where raw_counts measure the absolute number of reads and scaled_counts weights the read counts with the gene length. base.coverage measure the average base coverage of the gene.
> marker genes cutoff (-g). Every mOTU is composed of 10 marker genes and the read count of the mOTU is calculated as **the median of the read counts of the genes that are different from zero**. The parameter -g defines the minimum number of genes that have to be different from zero. The default value is 3 and possible values are between 1 and 10. With -g 1 the detection of one gene is enough to consider the mOTU as present in the sample (detecting low abundance species but also also false positives). On the other hand, with -g 6 only the mOTUs with 6 detected genes are counted, reducing the false positives.

```
>motus profile

Usage: motus profile [options]

Input options:
   -f FILE[,FILE]  input file(s) for reads in forward orientation, fastq formatted
   -r FILE[,FILE]  input file(s) for reads in reverse orientation, fastq formatted
   -s FILE[,FILE]  input file(s) for reads without mate, fastq formatted
   -n STR          sample name
   -i FILE[,FILE]  provide sam or bam input file(s)
   -m FILE         provide a mgc reads count file

Output options:
   -o FILE         output file name [stdout]
   -I FILE         save the result of bwa in bam format (intermediate step) [None]
   -M FILE         save the mgc reads count (intermediate step) [None]
   -e              profile only reference species (ref_mOTUs)
   -c              print result as counts instead of relative abundances
   -p              print NCBI id
   -u              print the full name of the species
   -q              print the full rank taxonomy
   -B              print result in BIOM format
   -C STR          print result in CAMI format (BioBoxes format 0.9.1)
                   Values: [precision, recall, parenthesis]
   -k STR          taxonomic level [mOTU]
                   Values: [kingdom, phylum, class, order, family, genus, mOTU]

Algorithm options:
   -g INT          number of marker genes cutoff: 1=higher recall, 10=higher precision [3]
   -l INT          min. length of alignment for the reads (number of nucleotides) [75]
   -t INT          number of threads [1]
   -v INT          verbose level: 1=error, 2=warning, 3=message, 4+=debugging [3]
   -y STR          type of read counts [insert.scaled_counts]
                   Values: [base.coverage, insert.raw_counts, insert.scaled_counts]
```

### 2. Merge profiles

```
>motus merge

Usage: motus merge [options]

Input options:
   -i STR[,STR]  list of files (comma separated)
   -d DIR        merge all the files in the directory

Output options:
   -o FILE       output file name [stdout]
   -B            print result in BIOM format

Algorithm options:
   -v INT        verbose level: 1=error, 2=warning, 3=message, 4+=debugging [3]
```

### 3. Output

- ヘッダー１：スクリプトとデータベースのバージョン情報（Ref-mOTUsはNCBIなどのレファレンスゲノムからのもの、meta-mOTUsは環境資料由来の多くは種名がついていないもの）
- ヘッダー２：コマンド
- ヘッダー３：列名
- unassigned：ファイルの最後にある。マップされなかったもの。定量はできない（特定の１種の値ではない）ので、解析時は使わないことを推奨。ただしrelative abundanceはunassignedを入れた数値をそのまま使うことを推奨（マップされた種を過大評価してしまうので）。

見つかるmOTUカウント数の目安（-cオプションで得られる値、human fecal sample）
|Total number of reads (million)|Median mOTUs count|
|---|---|
|5|600|
|8|900|
|15|1,900|
|25|3,300|
|35|5,500|
|50|8,800|
|100|13,000|

## SNV calling

### Map_snv
```
>motus map_snv

Usage: motus map_snv [options]

Input options:
   -f FILE[,FILE]  input file(s) for reads in forward orientation, fastq formatted
   -r FILE[,FILE]  input file(s) for reads in reverse orientation, fastq formatted
   -s FILE[,FILE]  input file(s) for reads without mate, fastq formatted

Output options:
   -o FILE         output bam file name [stdout]

Algorithm options:
   -l INT          min. length of alignment for the reads (number of nucleotides) [75]
   -t INT          number of threads [1]
   -v INT          verbose level: 1=error, 2=warning, 3=message, 4+=debugging [3]

```


### Snv_call

```
motus snv_call

Usage: motus snv_call -d Directory -o Directory [options]

Input options:
   -d  DIR     Call metaSNV on all bam files in the directory. [Mandatory]
   -fb FLOAT   Coverage breadth: minimal horizontal genome coverage percentage per sample per species. Default=80.0
   -fd FLOAT   Coverage depth: minimal average vertical genome coverage per sample per species. Default=5.0
   -fm INT     Minimum number of samples per species. Default=2
   -fp FLOAT   FILTERING STEP II: Required proportion of informative samples (coverage non-zero) per position. Default=0.50
   -fc FLOAT   FILTERING STEP II: Minimum coverage per position per sample per species. Default=5.0
   -t  INT     Number of threads. Default=1

Output options:
   -o  DIR     Output directory. Will fail if already exists. [Mandatory]
   -K          Keep all the directories produced by metaSNV. Default is to remove cov, distances, filtered, snpCaller

Algorithm options:
   -v INT      Verbose level: 1=error, 2=warning, 3=message, 4+=debugging. Default=3
```