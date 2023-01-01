# MicrobeCensus

- [MicrobeCensus wiki](https://github.com/snayfach/MicrobeCensus)

インストール
```
conda create -n MicrobeCensus -c bioconda microbecensus
```

テスト
```
cd ~/miniconda3/envs/MicrobeCensus/test
python test_microbe_census.py
```

ヘルプ
```
>run_microbe_census.py -h
usage: run_microbe_census.py [-h] [-v] [-V] [-r RAPSEARCH] [-n NREADS]
                             [-t THREADS] [-e]
                             [-l {50,60,70,80,90,100,110,120,130,140,150,175,200,225,250,300,350,400,450,500}]
                             [-q MIN_QUALITY] [-m MEAN_QUALITY] [-d]
                             [-u MAX_UNKNOWN]
                             SEQFILES OUTFILE

Estimate average genome size from metagenomic data.

positional arguments:
  SEQFILES              path to input metagenome(s); for paired-end
                        metagenomes use commas to specify each file (ex:
                        read_1.fq.gz,read_2.fq.gz); can be FASTQ/FASTA; can be
                        gzip (.gz) or bzip (.bz2) compressed
  OUTFILE               path to output file containing results

optional arguments:
  -h, --help            show this help message and exit
  -v                    print program's progress to stdout (default = False)
  -V, --version         show program's version number and exit
  -r RAPSEARCH          path to external RAPsearch2 v2.15 binary; useful if
                        precompiled RAPsearch2 v2.15 binary included with
                        MicrobeCensus does not work on your system

Pipeline throughput (optional):
  -n NREADS             number of reads to sample from SEQFILES and use for
                        average genome size estimation. to use all reads set
                        to 100000000. (default = 2000000)
  -t THREADS            number of threads to use for database search (default
                        = 1)
  -e                    quit after average genome size is obtained and do not
                        estimate the number of genome equivalents in SEQFILES.
                        useful in combination with -n for quick tests (default
                        = False)

Quality control (optional):
  -l {50,60,70,80,90,100,110,120,130,140,150,175,200,225,250,300,350,400,450,500}
                        all reads trimmed to this length; reads shorter than
                        this discarded (default = median read length)
  -q MIN_QUALITY        minimum base-level PHRED quality score (default = -5;
                        no filtering)
  -m MEAN_QUALITY       minimum read-level PHRED quality score (default = -5;
                        no filtering)
  -d                    filter duplicate reads (default = False)
  -u MAX_UNKNOWN        max percent of unknown bases per read (default = 100
                        percent; no filtering)
```