# nanopack

[wdecoster/nanopack](https://github.com/wdecoster/nanopack)

インストール
```
conda install -c bioconda nanoqc nanostat nanoplot nanofilt # など。condaだとバラで入れるみたい。
```


## nanoQC

htmlのグラフが出る
- リード長の分布
- 3'/5'末端からの位置 VS 塩基の頻度
- 3'/5'末端からの位置 VS 平均クオリティ

Example
```
nanoQC --outdir nanopack input.fastq.gz
```

Output
```
ls nanoqc
# nanoQC.html  NanoQC.log
```

Help
```
usage: nanoQC [-h] [-v] [-o OUTDIR] [--rna] [-l MINLEN] fastq

Investigate nucleotide composition and base quality.

positional arguments:
  fastq                 Reads data in fastq.gz format.

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         Print version and exit.
  -o OUTDIR, --outdir OUTDIR
                        Specify directory in which output has to be created.
  --rna                 Fastq is from direct RNA-seq and contains U nucleotides.
  -l MINLEN, --minlen MINLEN
                        Filters the reads on a minimal length of the given range. Also plots the
                        given length/2 of the begin and end of the reads.
```


## NanoStat

summaryのテキストファイル`NanoStats.txt`が作成される。標準出力に出るみたい。
下のNanoPlotでたくさん出てくる出力に含まれるので、NanoPlotをする場合はNanoStatはしなくていい。

Example
```
NanoStat -t 8 --fastq input.fastq.gz >nanopack/nanostats.txt
```

Help
```
usage: NanoStat [-h] [-v] [-o OUTDIR] [-p PREFIX] [-n NAME] [-t N] [--tsv] [--barcoded]
                [--readtype {1D,2D,1D2}] [--no_supplementary]
                (--fastq file [file ...] | --fasta file [file ...] | --summary file [file ...] | --bam file [file ...] | --ubam file [file ...] | --cram file [file ...] | --feather file [file ...])

Calculate statistics of long read sequencing dataset.

General options:
  -h, --help            show the help and exit
  -v, --version         Print version and exit.
  -o, --outdir OUTDIR   Specify directory for output, only in combination with -n.
  -p, --prefix PREFIX   Specify an optional prefix to be used for the output file.
  -n, --name NAME       Specify a filename/path for the output, stdout is the default.
  -t, --threads N       Set the allowed number of threads to be used by the script.
  --tsv                 Output the stats as a properly formatted TSV.

Input options.:
  --barcoded            Use if you want to split the summary file by barcode
  --readtype {1D,2D,1D2}
                        Which read type to extract information about from summary. Options are 1D, 2D,
                        1D2
  --no_supplementary    Use if you want to remove supplementary alignments

Input data sources, one of these is required.:
  --fastq file [file ...]
                        Data is in one or more (compressed) fastq file(s).
  --fasta file [file ...]
                        Data is in one or more (compressed) fasta file(s).
  --summary file [file ...]
                        Data is in one or more (compressed) summary file(s)generated by albacore or
                        guppy.
  --bam file [file ...]
                        Data is in one or more sorted bam file(s).
  --ubam file [file ...]
                        Data is in one or more unmapped bam file(s).
  --cram file [file ...]
                        Data is in one or more sorted cram file(s).
  --feather file [file ...]
                        Data is in one or more feather file(s).

EXAMPLES:
  NanoStat --fastq reads.fastq.gz --outdir statreports
  NanoStat --summary sequencing_summary1.txt sequencing_summary2.txtsequencing_summary3.txt --readtype 1D2
  NanoStat --bam alignment.bam alignment2.bam
```

## NanoPlot

ファイルがいっぱい出てくるのでこれ用のフォルダに入れよう。
NanoStatで出てくるsammaryも出力に含まれる。
prefixNanoPlot-report.htmlが各グラフをまとめたもの。

Example
```
NanoPlot -t 8 --loglength --fastq input.fastq.gz --outdir nanopack/nanoplot
```

Help
```
usage: NanoPlot [-h] [-v] [-t THREADS] [--verbose] [--store] [--raw] [--huge] [-o OUTDIR]
                [--no_static] [-p PREFIX] [--tsv_stats] [--info_in_report] [--maxlength N]
                [--minlength N] [--drop_outliers] [--downsample N] [--loglength] [--percentqual]
                [--alength] [--minqual N] [--runtime_until N] [--readtype {1D,2D,1D2}]
                [--barcoded] [--no_supplementary] [-c COLOR] [-cm COLORMAP]
                [-f [{png,jpg,jpeg,webp,svg,pdf,eps,json} [{png,jpg,jpeg,webp,svg,pdf,eps,json} ...]]]
                [--plots [{kde,hex,dot} [{kde,hex,dot} ...]]]
                [--legacy [{kde,dot,hex} [{kde,dot,hex} ...]]] [--listcolors] [--listcolormaps]
                [--no-N50] [--N50] [--title TITLE] [--font_scale FONT_SCALE] [--dpi DPI]
                [--hide_stats]
                (--fastq file [file ...] | --fasta file [file ...] | --fastq_rich file [file ...] | --fastq_minimal file [file ...] | --summary file [file ...] | --bam file [file ...] | --ubam file [file ...] | --cram file [file ...] | --pickle pickle | --feather file [file ...])

CREATES VARIOUS PLOTS FOR LONG READ SEQUENCING DATA.

General options:
  -h, --help            show the help and exit
  -v, --version         Print version and exit.
  -t, --threads THREADS
                        Set the allowed number of threads to be used by the script
  --verbose             Write log messages also to terminal.
  --store               Store the extracted data in a pickle file for future plotting.
  --raw                 Store the extracted data in tab separated file.
  --huge                Input data is one very large file.
  -o, --outdir OUTDIR   Specify directory in which output has to be created.
  --no_static           Do not make static (png) plots.
  -p, --prefix PREFIX   Specify an optional prefix to be used for the output files.
  --tsv_stats           Output the stats file as a properly formatted TSV.
  --info_in_report      Add NanoPlot run info in the report.

Options for filtering or transforming input prior to plotting:
  --maxlength N         Hide reads longer than length specified.
  --minlength N         Hide reads shorter than length specified.
  --drop_outliers       Drop outlier reads with extreme long length.
  --downsample N        Reduce dataset to N reads by random sampling.
  --loglength           Additionally show logarithmic scaling of lengths in plots.
  --percentqual         Use qualities as theoretical percent identities.
  --alength             Use aligned read lengths rather than sequenced length (bam mode)
  --minqual N           Drop reads with an average quality lower than specified.
  --runtime_until N     Only take the N first hours of a run
  --readtype {1D,2D,1D2}
                        Which read type to extract information about from summary. Options are 1D, 2D,
                        1D2
  --barcoded            Use if you want to split the summary file by barcode
  --no_supplementary    Use if you want to remove supplementary alignments

Options for customizing the plots created:
  -c, --color COLOR     Specify a valid matplotlib color for the plots
  -cm, --colormap COLORMAP
                        Specify a valid matplotlib colormap for the heatmap
  -f, --format [{png,jpg,jpeg,webp,svg,pdf,eps,json} [{png,jpg,jpeg,webp,svg,pdf,eps,json} ...]]
                        Specify the output format of the plots, which are in addition to the html files
  --plots [{kde,hex,dot} [{kde,hex,dot} ...]]
                        Specify which bivariate plots have to be made.
  --legacy [{kde,dot,hex} [{kde,dot,hex} ...]]
                        Specify which bivariate plots have to be made (legacy mode).
  --listcolors          List the colors which are available for plotting and exit.
  --listcolormaps       List the colors which are available for plotting and exit.
  --no-N50              Hide the N50 mark in the read length histogram
  --N50                 Show the N50 mark in the read length histogram
  --title TITLE         Add a title to all plots, requires quoting if using spaces
  --font_scale FONT_SCALE
                        Scale the font of the plots by a factor
  --dpi DPI             Set the dpi for saving images
  --hide_stats          Not adding Pearson R stats in some bivariate plots

Input data sources, one of these is required.:
  --fastq file [file ...]
                        Data is in one or more default fastq file(s).
  --fasta file [file ...]
                        Data is in one or more fasta file(s).
  --fastq_rich file [file ...]
                        Data is in one or more fastq file(s) generated by albacore, MinKNOW or guppy
                        with additional information concerning channel and time.
  --fastq_minimal file [file ...]
                        Data is in one or more fastq file(s) generated by albacore, MinKNOW or guppy
                        with additional information concerning channel and time. Is extracted swiftly
                        without elaborate checks.
  --summary file [file ...]
                        Data is in one or more summary file(s) generated by albacore or guppy.
  --bam file [file ...]
                        Data is in one or more sorted bam file(s).
  --ubam file [file ...]
                        Data is in one or more unmapped bam file(s).
  --cram file [file ...]
                        Data is in one or more sorted cram file(s).
  --pickle pickle       Data is a pickle file stored earlier.
  --feather file [file ...]
                        Data is in one or more feather file(s).

EXAMPLES:
    NanoPlot --summary sequencing_summary.txt --loglength -o summary-plots-log-transformed
    NanoPlot -t 2 --fastq reads1.fastq.gz reads2.fastq.gz --maxlength 40000 --plots hex dot
    NanoPlot --color yellow --bam alignment1.bam alignment2.bam alignment3.bam --downsample 10000
```

## NanoFilt

Example

```
gunzip -c reads.fastq.gz | NanoFilt -q 12 --headcrop 75 | gzip > trimmed-reads.fastq.gz
```

Help

```
usage: NanoFilt [-h] [-v] [--logfile LOGFILE] [-l LENGTH] [--maxlength MAXLENGTH] [-q QUALITY]
                [--minGC MINGC] [--maxGC MAXGC] [--headcrop HEADCROP] [--tailcrop TAILCROP]
                [-s SUMMARY] [--readtype {1D,2D,1D2}]
                [input]

Perform quality and/or length and/or GC filtering of (long read) fastq data.           Reads on stdin.

General options:
  -h, --help            show the help and exit
  -v, --version         Print version and exit.
  --logfile LOGFILE     Specify the path and filename for the log file.
  input                 input, uncompressed fastq file

Options for filtering reads on.:
  -l LENGTH, --length LENGTH
                        Filter on a minimum read length
  --maxlength MAXLENGTH
                        Filter on a maximum read length
  -q QUALITY, --quality QUALITY
                        Filter on a minimum average read quality score
  --minGC MINGC         Sequences must have GC content >= to this. Float between 0.0 and 1.0. Ignored if
                        using summary file.
  --maxGC MAXGC         Sequences must have GC content <= to this. Float between 0.0 and 1.0. Ignored if
                        using summary file.

Options for trimming reads.:
  --headcrop HEADCROP   Trim n nucleotides from start of read
  --tailcrop TAILCROP   Trim n nucleotides from end of read

Input options.:
  -s SUMMARY, --summary SUMMARY
                        Use albacore or guppy summary file for quality scores
  --readtype {1D,2D,1D2}
                        Which read type to extract information about from summary. Options are 1D, 2D or
                        1D2

EXAMPLES:
  gunzip -c reads.fastq.gz | NanoFilt -q 10 -l 500 --headcrop 50 | minimap2 genome.fa - | samtools sort -O BAM -@24 -o alignment.bam -
  gunzip -c reads.fastq.gz | NanoFilt -q 12 --headcrop 75 | gzip > trimmed-reads.fastq.gz
  gunzip -c reads.fastq.gz | NanoFilt -q 10 | gzip > highQuality-reads.fastq.gz
```

