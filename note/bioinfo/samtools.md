# Samtools

- [Samtools](http://www.htslib.org/)
- [Samtools - Manuals - samtools](http://www.htslib.org/doc/samtools.html)
- [Decoding SAM flags](https://broadinstitute.github.io/picard/explain-flags.html) SAMのフラグ

```sh
conda install -c bioconda samtools
```
samtools 1.13

samtools sort
```
# samファイルをソートしてbamファイルで出力
samtools sort -@ 4 -O BAM -o ${sample}.sort.bam ${sample}.sam 

>samtools sort
Usage: samtools sort [options...] [in.bam]
Options:
  -l INT     Set compression level, from 0 (uncompressed) to 9 (best)
  -u         Output uncompressed data (equivalent to -l 0)
  -m INT     Set maximum memory per thread; suffix K/M/G recognized [768M]
  -M         Use minimiser for clustering unaligned/unplaced reads
  -K INT     Kmer size to use for minimiser [20]
  -n         Sort by read name (not compatible with samtools index command)
  -t TAG     Sort by value of TAG. Uses position as secondary index (or read name if -n is set)
  -o FILE    Write final output to FILE rather than standard output
  -T PREFIX  Write temporary files to PREFIX.nnnn.bam
  --no-PG    do not add a PG line
      --input-fmt-option OPT[=VAL]
               Specify a single input file format option in the form
               of OPTION or OPTION=VALUE
  -O, --output-fmt FORMAT[,OPT[=VAL]]...
               Specify output format (SAM, BAM, CRAM)
      --output-fmt-option OPT[=VAL]
               Specify a single output file format option in the form
               of OPTION or OPTION=VALUE
      --reference FILE
               Reference sequence FASTA FILE [null]
  -@, --threads INT
               Number of additional threads to use [0]
      --write-index
               Automatically index the output files [off]
      --verbosity INT
               Set level of verbosity
```


### samtools view

[Samtools - samtools view](http://www.htslib.org/doc/samtools-view.html)

ヘルプ

```
Usage: samtools view [options] <in.bam>|<in.sam>|<in.cram> [region ...]

Output options:
  -b, --bam                  Output BAM
  -C, --cram                 Output CRAM (requires -T)
  -1, --fast                 Use fast BAM compression (implies --bam)
  -u, --uncompressed         Uncompressed BAM output (implies --bam)
  -h, --with-header          Include header in SAM output
  -H, --header-only          Print SAM header only (no alignments)
      --no-header            Print SAM alignment records only [default]
  -c, --count                Print only the count of matching records
  -o, --output FILE          Write output to FILE [standard output]
  -U, --unoutput FILE, --output-unselected FILE
                             Output reads not selected by filters to FILE
Input options:
  -t, --fai-reference FILE   FILE listing reference names and lengths
  -M, --use-index            Use index and multi-region iterator for regions
      --region[s]-file FILE  Use index to include only reads overlapping FILE
  -X, --customized-index     Expect extra index file argument after <in.bam>

Filtering options (Only include in output reads that...):
  -L, --target[s]-file FILE  ...overlap (BED) regions in FILE
  -r, --read-group STR       ...are in read group STR
  -R, --read-group-file FILE ...are in a read group listed in FILE
  -N, --qname-file FILE      ...whose read name is listed in FILE
  -d, --tag STR1[:STR2]      ...have a tag STR1 (with associated value STR2)
  -D, --tag-file STR:FILE    ...have a tag STR whose value is listed in FILE
  -q, --min-MQ INT           ...have mapping quality >= INT
  -l, --library STR          ...are in library STR
  -m, --min-qlen INT         ...cover >= INT query bases (as measured via CIGAR)
  -e, --expr STR             ...match the filter expression STR
  -f, --require-flags FLAG   ...have all of the FLAGs present
  -F, --excl[ude]-flags FLAG ...have none of the FLAGs present
  -G FLAG                    EXCLUDE reads with all of the FLAGs present
      --subsample FLOAT      Keep only FLOAT fraction of templates/read pairs
      --subsample-seed INT   Influence WHICH reads are kept in subsampling [0]
  -s INT.FRAC                Same as --subsample 0.FRAC --subsample-seed INT

Processing options:
      --add-flags FLAG       Add FLAGs to reads
      --remove-flags FLAG    Remove FLAGs from reads
  -x, --remove-tag STR       Strip tag STR from reads (option may be repeated)
  -B, --remove-B             Collapse the backward CIGAR operation

General options:
  -?, --help   Print long help, including note about region specification
  -S           Ignored (input format is auto-detected)
      --no-PG  Do not add a PG line
      --input-fmt-option OPT[=VAL]
               Specify a single input file format option in the form
               of OPTION or OPTION=VALUE
  -O, --output-fmt FORMAT[,OPT[=VAL]]...
               Specify output format (SAM, BAM, CRAM)
      --output-fmt-option OPT[=VAL]
               Specify a single output file format option in the form
               of OPTION or OPTION=VALUE
  -T, --reference FILE
               Reference sequence FASTA FILE [null]
  -@, --threads INT
               Number of additional threads to use [0]
      --write-index
               Automatically index the output files [off]
      --verbosity INT
               Set level of verbosity
```

R1もR2もマップされなかったリードのfastq

```
samtools view -h -f 12 ${sample}.sam \
| samtools fastq \
  -1 ${sample}_unmapped_R1.fastq \
  -2 ${sample}_unmapped_R2.fastq
```

R1とR2のいずれかがマップされたリードのfastq

```
samtools view -h -F 256 ${sample}.sam \
 | samtools view -h -G 12 \
 | samtools fastq \
  -1 ${sample}_mapped_R1.fastq \
  -2 ${sample}_mapped_R2.fastq
```


### 未整理

```
# secondary/supplementary alignmentを捨てる（１リード１結果になる）
samtools view -h -F 2048 sample.sam | samtools view -h -F 256 > sample_primary.sam

＃リード数をカウント
samtools view -c sample.sam
```

samtools sort
```
# 名前でsortしてbamにする
samtools sort -O BAM -n SAMPLE.sam -o SAMPLE.sort.bam
```
