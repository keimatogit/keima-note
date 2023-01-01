# spades


[github](https://github.com/ablab/spades#sec1)
[github - meta](https://github.com/ablab/spades#meta)
[metaSPAdes](http://cab.spbu.ru/software/meta-spades/)

```
conda install -c bioconda spades
```


```
SPAdes genome assembler v3.15.2

Usage: spades.py [options] -o <output_dir>

Basic options:
  -o <output_dir>             directory to store all the resulting files (required)
  --isolate                   this flag is highly recommended for high-coverage isolate and multi-cell data
  --sc                        this flag is required for MDA (single-cell) data
  --meta                      this flag is required for metagenomic data
  --bio                       this flag is required for biosyntheticSPAdes mode
  --corona                    this flag is required for coronaSPAdes mode
  --rna                       this flag is required for RNA-Seq data
  --plasmid                   runs plasmidSPAdes pipeline for plasmid detection
  --metaviral                 runs metaviralSPAdes pipeline for virus detection
  --metaplasmid               runs metaplasmidSPAdes pipeline for plasmid detection in metagenomic datasets (equivalent for --meta --plasmid)
  --rnaviral                  this flag enables virus assembly module from RNA-Seq data
  --iontorrent                this flag is required for IonTorrent data
  --test                      runs SPAdes on toy dataset
  -h, --help                  prints this usage message
  -v, --version               prints version

Input data:
  --12 <filename>             file with interlaced forward and reverse paired-end reads
  -1 <filename>               file with forward paired-end reads
  -2 <filename>               file with reverse paired-end reads
  -s <filename>               file with unpaired reads
  --merged <filename>         file with merged forward and reverse paired-end reads
  --pe-12 <#> <filename>      file with interlaced reads for paired-end library number <#>.
                              Older deprecated syntax is -pe<#>-12 <filename>
  --pe-1 <#> <filename>       file with forward reads for paired-end library number <#>.
                              Older deprecated syntax is -pe<#>-1 <filename>
  --pe-2 <#> <filename>       file with reverse reads for paired-end library number <#>.
                              Older deprecated syntax is -pe<#>-2 <filename>
  --pe-s <#> <filename>       file with unpaired reads for paired-end library number <#>.
                              Older deprecated syntax is -pe<#>-s <filename>
  --pe-m <#> <filename>       file with merged reads for paired-end library number <#>.
                              Older deprecated syntax is -pe<#>-m <filename>
  --pe-or <#> <or>            orientation of reads for paired-end library number <#> 
                              (<or> = fr, rf, ff).
                              Older deprecated syntax is -pe<#>-<or>
  --s <#> <filename>          file with unpaired reads for single reads library number <#>.
                              Older deprecated syntax is --s<#> <filename>
  --mp-12 <#> <filename>      file with interlaced reads for mate-pair library number <#>.
                              Older deprecated syntax is -mp<#>-12 <filename>
  --mp-1 <#> <filename>       file with forward reads for mate-pair library number <#>.
                              Older deprecated syntax is -mp<#>-1 <filename>
  --mp-2 <#> <filename>       file with reverse reads for mate-pair library number <#>.
                              Older deprecated syntax is -mp<#>-2 <filename>
  --mp-s <#> <filename>       file with unpaired reads for mate-pair library number <#>.
                              Older deprecated syntax is -mp<#>-s <filename>
  --mp-or <#> <or>            orientation of reads for mate-pair library number <#> 
                              (<or> = fr, rf, ff).
                              Older deprecated syntax is -mp<#>-<or>
  --hqmp-12 <#> <filename>    file with interlaced reads for high-quality mate-pair library number <#>.
                              Older deprecated syntax is -hqmp<#>-12 <filename>
  --hqmp-1 <#> <filename>     file with forward reads for high-quality mate-pair library number <#>.
                              Older deprecated syntax is -hqmp<#>-1 <filename>
  --hqmp-2 <#> <filename>     file with reverse reads for high-quality mate-pair library number <#>.
                              Older deprecated syntax is -hqmp<#>-2 <filename>
  --hqmp-s <#> <filename>     file with unpaired reads for high-quality mate-pair library number <#>.
                              Older deprecated syntax is -hqmp<#>-s <filename>
  --hqmp-or <#> <or>          orientation of reads for high-quality mate-pair library number <#> 
                              (<or> = fr, rf, ff).
                              Older deprecated syntax is -hqmp<#>-<or>
  --sanger <filename>         file with Sanger reads
  --pacbio <filename>         file with PacBio reads
  --nanopore <filename>       file with Nanopore reads
  --trusted-contigs <filename>
                              file with trusted contigs
  --untrusted-contigs <filename>
                              file with untrusted contigs

Pipeline options:
  --only-error-correction     runs only read error correction (without assembling)
  --only-assembler            runs only assembling (without read error correction)
  --careful                   tries to reduce number of mismatches and short indels
  --checkpoints <last or all>
                              save intermediate check-points ('last', 'all')
  --continue                  continue run from the last available check-point (only -o should be specified)
  --restart-from <cp>         restart run with updated options and from the specified check-point
                              ('ec', 'as', 'k<int>', 'mc', 'last')
  --disable-gzip-output       forces error correction not to compress the corrected reads
  --disable-rr                disables repeat resolution stage of assembling

Advanced options:
  --dataset <filename>        file with dataset description in YAML format
  -t <int>, --threads <int>   number of threads. [default: 16]
  -m <int>, --memory <int>    RAM limit for SPAdes in Gb (terminates if exceeded). [default: 250]
  --tmp-dir <dirname>         directory for temporary files. [default: <output_dir>/tmp]
  -k <int> [<int> ...]        list of k-mer sizes (must be odd and less than 128)
                              [default: 'auto']
  --cov-cutoff <float>        coverage cutoff value (a positive float number, or 'auto', or 'off')
                              [default: 'off']
  --phred-offset <33 or 64>   PHRED quality offset in the input reads (33 or 64),
                              [default: auto-detect]
  --custom-hmms <dirname>     directory with custom hmms that replace default ones,
                              [default: None]

```


## output

`-o`で指定したoutdirの中身：
- corrected/
- contigs.fasta
- scaffolds.fasta
- contigs.paths
- scaffolds.paths
- assembly_graph.fastg
- assembly_graph_with_scaffolds.gfa

[Difference Between Contig and Scaffold](https://www.differencebetween.com/difference-between-contig-and-scaffold/)
gene catalogue用にはcontigかな？

file             format  type  num_seqs      sum_len  min_len  avg_len  max_len
contigs.fasta    FASTA   DNA    401,266  332,415,778       56    828.4  695,975
scaffolds.fasta  FASTA   DNA    400,139  332,592,370       56    831.2  695,975


メタゲノムモードでは複数ファイルを指定してのco-assemblyはできない。co-assembleしたいときは１ファイルにまとめる？
[how to co-assemble multiple metagenome samples by using metaspades #656](https://github.com/ablab/spades/issues/656)

> --meta (same as metaspades.py) This flag is recommended when assembling metagenomic data sets (runs metaSPAdes, see paper for more details). Currently metaSPAdes supports only a single short-read library which has to be paired-end (we hope to remove this restriction soon). In addition, you can provide long reads (e.g. using --pacbio or --nanopore options), but hybrid assembly for metagenomes remains an experimental pipeline and optimal performance is not guaranteed. It does not support careful mode (mismatch correction is not available). In addition, you cannot specify coverage cutoff for metaSPAdes. Note that metaSPAdes might be very sensitive to presence of the technical sequences remaining in the data (most notably adapter readthroughs), please run quality control and pre-process your data accordingly.



220727_shotggun_F3618 spadesが使えなかった：
- spadesがペア認識しない。どうしても↓の警告が出る。
- unmappedを抜き出す前までならOKだったけど（qc後30万リードでテスト）
- samtools sort を挟んでもだめ。専用チャンネルでv3.15.5にしてもだめ。、samtool後はsort -n を挟んでもだめ。fastq_pairしてもだめ。samtoolsでとったリード名を参照してseqkit grepでfastqから抽出してもだめ（なぜ・・）
- --metaフラグをなくしてみたら、違う警告が出た。
- でも結果ファイルは普通に出てるっぽいけども・・
- 山野さんに聞いたら、samから抜き出すと参照配列にアライメント後だから向きが変わってるからかもということで、samtoolsでとったリード名を参照してseqkit grepでもとのfastq（qc後）から抜き出したけどダメだった。そういえば抜き出してるのはunmappなので向きは変わってないのかも？数リード見たところ、samtoolsで抜き出し後も向きは変わってなかった。山野さん的にはsamtoolsは何かとイマイチらしい。

```
# コマンド
sbatch -p c6420 -c22 -J spades -o test_spades2.out --wrap="\
spades.py \
 --meta \
 --threads 22 \
 -1 V01_qc_hg19unmapped_grep_R1.fastq.paired.fq \
 -2 V01_qc_hg19unmapped_grep_R2.fastq.paired.fq \
 -o spades_res"


# 標準出力の警告部分（--metaありのとき）

======= SPAdes pipeline finished WITH WARNINGS!

=== Error correction and assembling warnings:
 * 0:00:06.480    15M / 422M  WARN    General                 (pair_info_count.cpp       : 358)   Unable to estimate insert size for paired library #0
 * 0:00:06.480    15M / 422M  WARN    General                 (pair_info_count.cpp       : 364)   None of paired reads aligned properly. Please, check orientation of your read pairs.
 * 0:00:06.481    15M / 422M  WARN    General                 (repeat_resolving.cpp      :  81)   Insert size was not estimated for any of the paired libraries, repeat resolution module will not run.
 * 0:00:07.768    12M / 542M  WARN    General                 (pair_info_count.cpp       : 187)   Single reads are not used in metagenomic mode
 
 
 # 標準出力の警告部分（--metaなしのとき）
 
 ======= SPAdes pipeline finished WITH WARNINGS!

=== Error correction and assembling warnings:
 * 0:00:07.042    10M / 158M  WARN    General                 (kmer_coverage_model.cpp   : 327)   Valley value was estimated improperly, reset to 9
 * 0:00:07.127     9M / 158M  WARN    General                 (simplification.cpp        : 500)   The determined erroneous connection coverage threshold may be determined improperly
 * 0:00:09.684    10M / 154M  WARN    General                 (kmer_coverage_model.cpp   : 327)   Valley value was estimated improperly, reset to 3
 * 0:00:09.685    10M / 154M  WARN    General                 (kmer_coverage_model.cpp   : 366)   Failed to determine erroneous kmer threshold. Threshold set to: 3
 * 0:00:05.270    10M / 102M  WARN    General                 (kmer_coverage_model.cpp   : 218)   Too many erroneous kmers, the estimates might be unreliable
 * 0:00:06.830    10M / 102M  WARN    General                 (kmer_coverage_model.cpp   : 327)   Valley value was estimated improperly, reset to 11
```

220807_shotgunで動いた。4時間半くらい。（--memory はデフォが250。これでメモリ不足になる場合も。）

```
sbatch -p c6420 -c22 -J spades -o slurm_spades_test_%j.out --wrap="\
spades.py \
 --threads 22 \
 --memory 250 \
 --meta \
 -1 FM9_R1.fastq.gz \
 -2 FM9_R2.fastq.gz \
 -o spades_test_res"
```
