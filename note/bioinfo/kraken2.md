# Kraken2, Bracken


[Metagenome analysis using the Kraken software suite](https://www.nature.com/articles/s41596-022-00738-y)

Kraken2
[DerrickWood/kraken2 - Manual](https://github.com/DerrickWood/kraken2/wiki/Manual)
[Kraken taxonomic sequence classification system](https://software.cqls.oregonstate.edu/updates/docs/kraken2/MANUAL.html)
[KrakenTools](https://github.com/jenniferlu717/KrakenTools)

Bracken
[Bracken - Bayesian Reestimation of Abundance with Kraken](https://ccb.jhu.edu/software/bracken/index.shtml?t=manual)
[jenniferlu717/Bracken](https://github.com/jenniferlu717/Bracken)
[kraken2/bracken setup and read-classification example](https://hackmd.io/@astrobiomike/kraken2-bracken-standard-build)


## インストール〜データベース作成

```sh
conda create -n kraken2 -c bioconda kraken2 bracken
```

```
# create standard Kraken 2 database
[srun -c 88] kraken2-build --standard [--threads 88] --db ~/db/kraken2/standard

# NCBI nt
srun -c 88 kraken2-build --threads 88 --download-library nt --db ~/db/kraken2/nt

# NCBI nr
srun -c 88 kraken2-build --threads 88 --protein --download-library nr --db ~/db/kraken2/nr
```

2022/03/17：kraken2-buildでエラー`rsync_from_ncbi.pl: unexpected FTP path (new server?) for https:// ... `が出たので、
githubから新しいrsync_from_ncbi.pl(https://github.com/DerrickWood/kraken2/blob/master/scripts/rsync_from_ncbi.pl)をコピペして、さらに`$full_path =~ s#^ftp://` => `$full_path =~ s#^https://`に修正したものを~/miniconda3/envs/kraken2/libexec/rsync_from_ncbi.plとして保存。パーミッション775もつけておくこと。

[Kraken2 でウイルス配列DBを取得するときのエラー対処](https://qiita.com/kohei-108/items/ce5fdf10c11d1e7ca15b)

```
rm ~/miniconda3/envs/kraken2/libexec/rsync_from_ncbi.pl
vi ~/miniconda3/envs/kraken2/libexec/rsync_from_ncbi.pl #githubのrsync_from_ncbi.plをコピペ、編集。
chmod 775 ~/miniconda3/envs/kraken2/libexec/rsync_from_ncbi.pl
```





## ヘルプ

kraken2 -h
```sh
Usage: kraken2 [options] <filename(s)>

Options:
  --db NAME               Name for Kraken 2 DB
                          (default: none)
  --threads NUM           Number of threads (default: 1)
  --quick                 Quick operation (use first hit or hits)
  --unclassified-out FILENAME
                          Print unclassified sequences to filename
  --classified-out FILENAME
                          Print classified sequences to filename
  --output FILENAME       Print output to filename (default: stdout); "-" will
                          suppress normal output
  --confidence FLOAT      Confidence score threshold (default: 0.0); must be
                          in [0, 1].
  --minimum-base-quality NUM
                          Minimum base quality used in classification (def: 0,
                          only effective with FASTQ input).
  --report FILENAME       Print a report with aggregrate counts/clade to file
  --use-mpa-style         With --report, format report output like Kraken 1's
                          kraken-mpa-report
  --report-zero-counts    With --report, report counts for ALL taxa, even if
                          counts are zero
  --report-minimizer-data With --report, report minimizer and distinct minimizer
                          count information in addition to normal Kraken report
  --memory-mapping        Avoids loading database into RAM
  --paired                The filenames provided have paired-end reads
  --use-names             Print scientific names instead of just taxids
  --gzip-compressed       Input files are compressed with gzip
  --bzip2-compressed      Input files are compressed with bzip2
  --minimum-hit-groups NUM
                          Minimum number of hit groups (overlapping k-mers
                          sharing the same minimizer) needed to make a call
                          (default: 2)
  --help                  Print this message
  --version               Print version information

If none of the *-compressed flags are specified, and the filename provided
is a regular file, automatic format detection is attempted.
```


## Example

kronaには--output、blackenには--reportのファイルが必要。

```bash
#!/bin/bash
#SBATCH -c 44

kraken2 \
 --threads 44 \
 --db  ~/db/kraken2/standard \
 --paired \
 --gzip-compressed \
 --output ${sample}_kraken.txt \
 --report ${sample}_kraken_report.txt \
 --use-mpa-style \
 ../0_fastq/${sample}_clean_R1.fastq.gz \
 ../0_fastq/${sample}_clean_R2.fastq.gz
```

## Output

### --output

> Each sequence classified by Kraken results in a single line of output. Output lines contain five tab-delimited fields; from left to right, they are:

1. "C"または"U": ClassifiedかUnclassifiedかを示す。
2. リード名。
3. taxonomy id (unclassifiedの場合は0)。
4. リード長(bp)。
5. A space-delimited list indicating the LCA mapping of each k-mer in the sequence. For example, "562:13 561:4 A:31 0:1 562:3" 

```
C       7b304bd4-e876-4ebf-9073-37862ca5d314    718255  917     718255:22 0:51 186803:8 718255:6 186803:3 718255:31 186803:1 718255:4 0:45 718255:18 0:39 718255:15 0:27 718255:22 0:1 718255:5 0:48 718255:10 0:5 718255:7 0:22 718255:13 1453429:2 0:101 718255:55 0:134 718255:5 0:58 718255:5 0:120
C       707ab901-a88f-43ba-9d49-868e752a0da4    821     744     0:9 171549:8 0:33 815:5 0:31 815:15 171549:5 815:1 171549:20 821:6 0:31 821:3 0:38 821:17 171549:1 821:10 171549:29 0:33 171549:8 815:3 171549:5 815:5 0:34 815:33 171549:6 815:22 0:39 171549:24 815:7 171549:2 815:5 171549:5 815:13 171549:9 815:33 171549:13 0:26 821:11 0:75 357276:37
C       cbe98012-5329-4ab5-8dcb-cde8a3c31796    562     3527    0:51 1:38 0:32 1:4 59201:5 0:73 1:14 0:100 543:1 562:1 1:11 2:3 1:7 2:6 0:54 1:26 0:15 1:3 0:8 1:12 0:88 1:119 0:40 1:8 0:29 1:14 0:48 1:109 0:42 1:78 0:23 1:1 0:7 1:62 562:1 1:10 562:5 1:8 562:5 1:50 562:7 1:5 0:42 1:23 562:15 0:5 562:2 561:5 1243595:9 91347:6 543:5 562:3 1:52 562:3 0:35 562:23 1:128 0:30 1:103 562:5 0:27 562:7 1:113 0:26 1:74 0:53 1:216 0:57 1:47 0:72 1:30 0:32 1:110 0:35 1:23 0:31 1:60 0:36 1:6 0:6 1:44 0:59 1:150 0:5 1:3 0:11 1:7 0:1 1:7 0:32 1:58 562:2 0:37 1:250 0:49
C       67899b1f-f042-4c73-a9b1-1e6807ef6df1    2       932     0:49 2780074:1 0:803 39485:1 0:44
U       6ec7366d-f8d9-471c-892e-3d6b68b855bf    0       877     0:843
```

### --report


> The output of kraken-report is tab-delimited, with one line per taxon. The fields of the output, from left-to-right, are as follows:

1. Percentage of reads covered by the clade rooted at this taxon
2. Number of reads covered by the clade rooted at this taxon
3. Number of reads assigned directly to this taxon
4. A rank code, indicating (U)nclassified, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. All other ranks are simply '-'.
5. NCBI taxonomy ID
6. indented scientific name（インデントはspace2つ単位）

```
 18.03	374508	374508	U	0	unclassified
 81.97	1702359	1075	R	1	root
 81.87	1700265	8596	R1	131567	  cellular organisms
 80.60	1674016	41629	D	2	    Bacteria
 37.05	769483	12	D1	1783270	      FCB group
 37.05	769380	23	D2	68336	        Bacteroidetes/Chlorobi group
 37.04	769206	1607	P	976	          Bacteroidetes
 36.58	759689	8	C	200643	            Bacteroidia
 36.57	759563	21022	O	171549	              Bacteroidales
 28.21	585940	26507	F	815	                Bacteroidaceae
 19.42	403275	27360	G	816	                  Bacteroides
  7.83	162715	116882	S	820	                    Bacteroides uniformis
  2.21	45833	45833	S1	997890	                      Bacteroides uniformis CL03T12C37
  3.62	75273	68938	S	371601	                    Bacteroides xylanisolvens
  0.31	6335	6335	S1	657309	                      Bacteroides xylanisolvens XB1A
  2.63	54602	54196	S	818	                    Bacteroides thetaiotaomicron
  0.02	406	406	S1	226186	                      Bacteroides thetaiotaomicron VPI-5482
```

## kronaでhtmlのグラフにする（複数ファイル可）

```
ktImportTaxonomy \
 -q 2 -t 3 ${sample}_kraken.txt -o ${sample}_kraken.html
```

## blacken

### blacken用のDB作成

`bracken`を実行する際に`-r`で指定するREAD_LENと同じREAD_LENを指定して、blacken用のDBを作成しておく必要がある。22スレッドで1時間半くらいかかった。

```
Usage: bracken_build -k KMER_LEN -l READ_LEN -d MY_DB -x K_INSTALLATION -t THREADS
  KMER_LEN       kmer length used to build the kraken database (default: 35)
  THREADS        the number of threads to use when running kraken classification and the bracken scripts
  READ_LEN       read length to get all classifications for (default: 100)
  MY_DB          location of Kraken database
  K_INSTALLATION location of the installed kraken/kraken-build scripts (default assumes scripts can be run from the user path)

**Note that this script will try to use kraken2 as default. If kraken2 is not installed, kraken will be used instead
```

Example

```
bracken-build \
  -t 22 \
  -l 120 \
  -d ~/db/kraken2/standard
```

例えば`-l 120`とした場合、database120mers.kraken、database120mers.kmer_distribという２ファイルが`-d`で指定したデータベースと同じフォルダ内に作成される。


### brackenの実行

```
Usage: bracken -d MY_DB -i INPUT -o OUTPUT -w OUTREPORT -r READ_LEN -l LEVEL -t THRESHOLD
  MY_DB          location of Kraken database
  INPUT          Kraken REPORT file to use for abundance estimation
  OUTPUT         file name for Bracken default output
  OUTREPORT      New Kraken REPORT output file with Bracken read estimates
  READ_LEN       read length to get all classifications for (default: 100)
  LEVEL          level to estimate abundance at [options: D,P,C,O,F,G,S,S1,etc] (default: S)
  THRESHOLD      number of reads required PRIOR to abundance estimation to perform reestimation (default: 0)
```

```
bracken \
 -d ~/db/kraken2/standard \
 -t 22 \
 -l G \
 -r 100 \
 -i ${sample}_kraken_report.txt \
 -o ${sample}_kraken.blacken
```


## KrakenTools

[https://github.com/jenniferlu717/KrakenTools](https://github.com/jenniferlu717/KrakenTools)

### extract_kraken_reads.py

```
> extract_kraken_reads.py -h

usage: extract_kraken_reads.py [-h] -k KRAKEN_FILE -s SEQ_FILE1 [-s2 SEQ_FILE2] -t TAXID
                               [TAXID ...] -o OUTPUT_FILE [-o2 OUTPUT_FILE2] [--append]
                               [--noappend] [--max MAX_READS] [-r REPORT_FILE] [--include-parents]
                               [--include-children] [--exclude] [--fastq-output]

optional arguments:
  -h, --help            show this help message and exit
  -k KRAKEN_FILE        Kraken output file to parse
  -s SEQ_FILE1, -s1 SEQ_FILE1, -1 SEQ_FILE1, -U SEQ_FILE1
                        FASTA/FASTQ File containing the raw sequence letters.
  -s2 SEQ_FILE2, -2 SEQ_FILE2
                        2nd FASTA/FASTQ File containing the raw sequence letters (paired).
  -t TAXID [TAXID ...], --taxid TAXID [TAXID ...]
                        Taxonomy ID[s] of reads to extract (space-delimited)
  -o OUTPUT_FILE, --output OUTPUT_FILE
                        Output FASTA/Q file containing the reads and sample IDs
  -o2 OUTPUT_FILE2, --output2 OUTPUT_FILE2
                        Output FASTA/Q file containig the second pair of reads [required for
                        paired input]
  --append              Append the sequences to the end of the output FASTA file specified.
  --noappend            Create a new FASTA file containing sample sequences and IDs (rewrite if
                        existing) [default].
  --max MAX_READS       Maximum number of reads to save [default: 100,000,000]
  -r REPORT_FILE, --report REPORT_FILE
                        Kraken report file. [required only if --include-parents/children is
                        specified]
  --include-parents     Include reads classified at parent levels of the specified taxids
  --include-children    Include reads classified more specifically than the specified taxids
  --exclude             Instead of finding reads matching specified taxids, finds all reads NOT
                        matching specified taxids
  --fastq-output        Print output FASTQ reads [requires input FASTQ, default: output is FASTA]
```


Example

```
extract_kraken_reads.py \
 --taxid 9606 \
 --exclude \
 --include-children \
 --fastq-output \
 -k  FM9_qc_kraken.txt \
 -r  FM9_qc_kraken.report \
 -s  FM9_qc_R1.fastq.gz \
 -s2 FM9_qc_R2.fastq.gz \
 -o  FM9_qc_noHuman_R1.fastq \
 -o2 FM9_qc_noHuman_R2.fastq

```# Kraken2, Bracken


[Metagenome analysis using the Kraken software suite](https://www.nature.com/articles/s41596-022-00738-y)

Kraken2
[DerrickWood/kraken2 - Manual](https://github.com/DerrickWood/kraken2/wiki/Manual)
[Kraken taxonomic sequence classification system](https://software.cqls.oregonstate.edu/updates/docs/kraken2/MANUAL.html)
[KrakenTools](https://github.com/jenniferlu717/KrakenTools)

Bracken
[Bracken - Bayesian Reestimation of Abundance with Kraken](https://ccb.jhu.edu/software/bracken/index.shtml?t=manual)
[jenniferlu717/Bracken](https://github.com/jenniferlu717/Bracken)
[kraken2/bracken setup and read-classification example](https://hackmd.io/@astrobiomike/kraken2-bracken-standard-build)


## インストール〜データベース作成

```sh
conda create -n kraken2 -c bioconda kraken2 bracken
```

```
# create standard Kraken 2 database
[srun -c 88] kraken2-build --standard [--threads 88] --db ~/db/kraken2/standard

# NCBI nt
srun -c 88 kraken2-build --threads 88 --download-library nt --db ~/db/kraken2/nt

# NCBI nr
srun -c 88 kraken2-build --threads 88 --protein --download-library nr --db ~/db/kraken2/nr
```

2022/03/17：kraken2-buildでエラー`rsync_from_ncbi.pl: unexpected FTP path (new server?) for https:// ... `が出たので、
githubから新しいrsync_from_ncbi.pl(https://github.com/DerrickWood/kraken2/blob/master/scripts/rsync_from_ncbi.pl)をコピペして、さらに`$full_path =~ s#^ftp://` => `$full_path =~ s#^https://`に修正したものを~/miniconda3/envs/kraken2/libexec/rsync_from_ncbi.plとして保存。パーミッション775もつけておくこと。

[Kraken2 でウイルス配列DBを取得するときのエラー対処](https://qiita.com/kohei-108/items/ce5fdf10c11d1e7ca15b)

```
rm ~/miniconda3/envs/kraken2/libexec/rsync_from_ncbi.pl
vi ~/miniconda3/envs/kraken2/libexec/rsync_from_ncbi.pl #githubのrsync_from_ncbi.plをコピペ、編集。
chmod 775 ~/miniconda3/envs/kraken2/libexec/rsync_from_ncbi.pl
```





## ヘルプ

kraken2 -h
```sh
Usage: kraken2 [options] <filename(s)>

Options:
  --db NAME               Name for Kraken 2 DB
                          (default: none)
  --threads NUM           Number of threads (default: 1)
  --quick                 Quick operation (use first hit or hits)
  --unclassified-out FILENAME
                          Print unclassified sequences to filename
  --classified-out FILENAME
                          Print classified sequences to filename
  --output FILENAME       Print output to filename (default: stdout); "-" will
                          suppress normal output
  --confidence FLOAT      Confidence score threshold (default: 0.0); must be
                          in [0, 1].
  --minimum-base-quality NUM
                          Minimum base quality used in classification (def: 0,
                          only effective with FASTQ input).
  --report FILENAME       Print a report with aggregrate counts/clade to file
  --use-mpa-style         With --report, format report output like Kraken 1's
                          kraken-mpa-report
  --report-zero-counts    With --report, report counts for ALL taxa, even if
                          counts are zero
  --report-minimizer-data With --report, report minimizer and distinct minimizer
                          count information in addition to normal Kraken report
  --memory-mapping        Avoids loading database into RAM
  --paired                The filenames provided have paired-end reads
  --use-names             Print scientific names instead of just taxids
  --gzip-compressed       Input files are compressed with gzip
  --bzip2-compressed      Input files are compressed with bzip2
  --minimum-hit-groups NUM
                          Minimum number of hit groups (overlapping k-mers
                          sharing the same minimizer) needed to make a call
                          (default: 2)
  --help                  Print this message
  --version               Print version information

If none of the *-compressed flags are specified, and the filename provided
is a regular file, automatic format detection is attempted.
```


## Example

kronaには--output、blackenには--reportのファイルが必要。

```bash
#!/bin/bash
#SBATCH -c 44

kraken2 \
 --threads 44 \
 --db  ~/db/kraken2/standard \
 --paired \
 --gzip-compressed \
 --output ${sample}_kraken.txt \
 --report ${sample}_kraken_report.txt \
 --use-mpa-style \
 ../0_fastq/${sample}_clean_R1.fastq.gz \
 ../0_fastq/${sample}_clean_R2.fastq.gz
```

## Output

### --output

> Each sequence classified by Kraken results in a single line of output. Output lines contain five tab-delimited fields; from left to right, they are:

1. "C"または"U": ClassifiedかUnclassifiedかを示す。
2. リード名。
3. taxonomy id (unclassifiedの場合は0)。
4. リード長(bp)。
5. A space-delimited list indicating the LCA mapping of each k-mer in the sequence. For example, "562:13 561:4 A:31 0:1 562:3" 

```
C       7b304bd4-e876-4ebf-9073-37862ca5d314    718255  917     718255:22 0:51 186803:8 718255:6 186803:3 718255:31 186803:1 718255:4 0:45 718255:18 0:39 718255:15 0:27 718255:22 0:1 718255:5 0:48 718255:10 0:5 718255:7 0:22 718255:13 1453429:2 0:101 718255:55 0:134 718255:5 0:58 718255:5 0:120
C       707ab901-a88f-43ba-9d49-868e752a0da4    821     744     0:9 171549:8 0:33 815:5 0:31 815:15 171549:5 815:1 171549:20 821:6 0:31 821:3 0:38 821:17 171549:1 821:10 171549:29 0:33 171549:8 815:3 171549:5 815:5 0:34 815:33 171549:6 815:22 0:39 171549:24 815:7 171549:2 815:5 171549:5 815:13 171549:9 815:33 171549:13 0:26 821:11 0:75 357276:37
C       cbe98012-5329-4ab5-8dcb-cde8a3c31796    562     3527    0:51 1:38 0:32 1:4 59201:5 0:73 1:14 0:100 543:1 562:1 1:11 2:3 1:7 2:6 0:54 1:26 0:15 1:3 0:8 1:12 0:88 1:119 0:40 1:8 0:29 1:14 0:48 1:109 0:42 1:78 0:23 1:1 0:7 1:62 562:1 1:10 562:5 1:8 562:5 1:50 562:7 1:5 0:42 1:23 562:15 0:5 562:2 561:5 1243595:9 91347:6 543:5 562:3 1:52 562:3 0:35 562:23 1:128 0:30 1:103 562:5 0:27 562:7 1:113 0:26 1:74 0:53 1:216 0:57 1:47 0:72 1:30 0:32 1:110 0:35 1:23 0:31 1:60 0:36 1:6 0:6 1:44 0:59 1:150 0:5 1:3 0:11 1:7 0:1 1:7 0:32 1:58 562:2 0:37 1:250 0:49
C       67899b1f-f042-4c73-a9b1-1e6807ef6df1    2       932     0:49 2780074:1 0:803 39485:1 0:44
U       6ec7366d-f8d9-471c-892e-3d6b68b855bf    0       877     0:843
```

### --report


> The output of kraken-report is tab-delimited, with one line per taxon. The fields of the output, from left-to-right, are as follows:

1. Percentage of reads covered by the clade rooted at this taxon
2. Number of reads covered by the clade rooted at this taxon
3. Number of reads assigned directly to this taxon
4. A rank code, indicating (U)nclassified, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. All other ranks are simply '-'.
5. NCBI taxonomy ID
6. indented scientific name（インデントはspace2つ単位）

```
 18.03	374508	374508	U	0	unclassified
 81.97	1702359	1075	R	1	root
 81.87	1700265	8596	R1	131567	  cellular organisms
 80.60	1674016	41629	D	2	    Bacteria
 37.05	769483	12	D1	1783270	      FCB group
 37.05	769380	23	D2	68336	        Bacteroidetes/Chlorobi group
 37.04	769206	1607	P	976	          Bacteroidetes
 36.58	759689	8	C	200643	            Bacteroidia
 36.57	759563	21022	O	171549	              Bacteroidales
 28.21	585940	26507	F	815	                Bacteroidaceae
 19.42	403275	27360	G	816	                  Bacteroides
  7.83	162715	116882	S	820	                    Bacteroides uniformis
  2.21	45833	45833	S1	997890	                      Bacteroides uniformis CL03T12C37
  3.62	75273	68938	S	371601	                    Bacteroides xylanisolvens
  0.31	6335	6335	S1	657309	                      Bacteroides xylanisolvens XB1A
  2.63	54602	54196	S	818	                    Bacteroides thetaiotaomicron
  0.02	406	406	S1	226186	                      Bacteroides thetaiotaomicron VPI-5482
```

## kronaでhtmlのグラフにする（複数ファイル可）

```
ktImportTaxonomy \
 -q 2 -t 3 ${sample}_kraken.txt -o ${sample}_kraken.html
```

## blacken

### blacken用のDB作成

`bracken`を実行する際に`-r`で指定するREAD_LENと同じREAD_LENを指定して、blacken用のDBを作成しておく必要がある。22スレッドで1時間半くらいかかった。

```
Usage: bracken_build -k KMER_LEN -l READ_LEN -d MY_DB -x K_INSTALLATION -t THREADS
  KMER_LEN       kmer length used to build the kraken database (default: 35)
  THREADS        the number of threads to use when running kraken classification and the bracken scripts
  READ_LEN       read length to get all classifications for (default: 100)
  MY_DB          location of Kraken database
  K_INSTALLATION location of the installed kraken/kraken-build scripts (default assumes scripts can be run from the user path)

**Note that this script will try to use kraken2 as default. If kraken2 is not installed, kraken will be used instead
```

Example

```
bracken-build \
  -t 22 \
  -l 120 \
  -d ~/db/kraken2/standard
```

例えば`-l 120`とした場合、database120mers.kraken、database120mers.kmer_distribという２ファイルが`-d`で指定したデータベースと同じフォルダ内に作成される。


### brackenの実行

```
Usage: bracken -d MY_DB -i INPUT -o OUTPUT -w OUTREPORT -r READ_LEN -l LEVEL -t THRESHOLD
  MY_DB          location of Kraken database
  INPUT          Kraken REPORT file to use for abundance estimation
  OUTPUT         file name for Bracken default output
  OUTREPORT      New Kraken REPORT output file with Bracken read estimates
  READ_LEN       read length to get all classifications for (default: 100)
  LEVEL          level to estimate abundance at [options: D,P,C,O,F,G,S,S1,etc] (default: S)
  THRESHOLD      number of reads required PRIOR to abundance estimation to perform reestimation (default: 0)
```

```
bracken \
 -d ~/db/kraken2/standard \
 -t 22 \
 -l G \
 -r 100 \
 -i ${sample}_kraken_report.txt \
 -o ${sample}_kraken.blacken
```


## KrakenTools

[https://github.com/jenniferlu717/KrakenTools](https://github.com/jenniferlu717/KrakenTools)

### extract_kraken_reads.py

```
> extract_kraken_reads.py -h

usage: extract_kraken_reads.py [-h] -k KRAKEN_FILE -s SEQ_FILE1 [-s2 SEQ_FILE2] -t TAXID
                               [TAXID ...] -o OUTPUT_FILE [-o2 OUTPUT_FILE2] [--append]
                               [--noappend] [--max MAX_READS] [-r REPORT_FILE] [--include-parents]
                               [--include-children] [--exclude] [--fastq-output]

optional arguments:
  -h, --help            show this help message and exit
  -k KRAKEN_FILE        Kraken output file to parse
  -s SEQ_FILE1, -s1 SEQ_FILE1, -1 SEQ_FILE1, -U SEQ_FILE1
                        FASTA/FASTQ File containing the raw sequence letters.
  -s2 SEQ_FILE2, -2 SEQ_FILE2
                        2nd FASTA/FASTQ File containing the raw sequence letters (paired).
  -t TAXID [TAXID ...], --taxid TAXID [TAXID ...]
                        Taxonomy ID[s] of reads to extract (space-delimited)
  -o OUTPUT_FILE, --output OUTPUT_FILE
                        Output FASTA/Q file containing the reads and sample IDs
  -o2 OUTPUT_FILE2, --output2 OUTPUT_FILE2
                        Output FASTA/Q file containig the second pair of reads [required for
                        paired input]
  --append              Append the sequences to the end of the output FASTA file specified.
  --noappend            Create a new FASTA file containing sample sequences and IDs (rewrite if
                        existing) [default].
  --max MAX_READS       Maximum number of reads to save [default: 100,000,000]
  -r REPORT_FILE, --report REPORT_FILE
                        Kraken report file. [required only if --include-parents/children is
                        specified]
  --include-parents     Include reads classified at parent levels of the specified taxids
  --include-children    Include reads classified more specifically than the specified taxids
  --exclude             Instead of finding reads matching specified taxids, finds all reads NOT
                        matching specified taxids
  --fastq-output        Print output FASTQ reads [requires input FASTQ, default: output is FASTA]
```


Example

```
extract_kraken_reads.py \
 --taxid 9606 \
 --exclude \
 --include-children \
 --fastq-output \
 -k  FM9_qc_kraken.txt \
 -r  FM9_qc_kraken.report \
 -s  FM9_qc_R1.fastq.gz \
 -s2 FM9_qc_R2.fastq.gz \
 -o  FM9_qc_noHuman_R1.fastq \
 -o2 FM9_qc_noHuman_R2.fastq

```