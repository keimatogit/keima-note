# tRNAscan

tRNAを探す

[github](https://github.com/UCSC-LoweLab/tRNAscan-SE)
[macでインフォマティクス- tRNAをゲノムやraw fastqから探す tRNAscan-SE 2.0](https://kazumaxneo.hatenablog.com/entry/2019/05/07/073000)

```
conda create -n trnascan -c bioconda trnascan-se
```

tRNAscan-SE -h

```
tRNAscan-SE 2.0.11 (Oct 2022)
Copyright (C) 2022 Patricia Chan and Todd Lowe
                   University of California Santa Cruz
Freely distributed under the GNU General Public License (GPLv3)


Usage: tRNAscan-SE [-options] <FASTA file(s)>

  Scan a sequence file for tRNAs 
   -- default: use Infernal & tRNA covariance models
      with eukaryotic sequences 
      (use 'Search Mode Options' below to scan other types of sequences)

Search Mode Options:

  -E                          : search for eukaryotic tRNAs (default)
  -B                          : search for bacterial tRNAs
  -A                          : search for archaeal tRNAs
  -M <model>                  : search for mitochondrial tRNAs
                                  options: mammal, vert
  -O                          : search for other organellar tRNAs
  -G                          : use general tRNA model (cytoslic tRNAs from all 3 domains included)
  --mt <model>                : use mito tRNA models for cytosolic/mito detemination
                                  (if not specified, only cytosolic isotype-specific model scan will be performed)
  -I                          : search using Infernal
                                  default use with -E, -B, -A, or -G; optional for -O
      --max                   : maximum sensitivity mode - search using Infernal without hmm filter (very slow)
  -L                          : search using the legacy method (tRNAscan, EufindtRNA, and COVE)
                                  use with -E, -B, -A or -G
  -C  --cove                  : search using COVE analysis only (legacy, extremely slow)
                                  default use with -O
  -H  --breakdown             : show breakdown of primary and secondary structure components to
                                  covariance model bit scores
  -D  --nopseudo              : disable pseudogene checking

Output options:

  -o  --output <file>         : save final results in <file>
  -f  --struct <file>         : save tRNA secondary structures to <file>
  -s  --isospecific <file>    : save results using isotype-specific models in <file>
  -m  --stats <file>          : save statistics summary for run in <file>
                                  (speed, # tRNAs found in each part of search, etc)
  -b  --bed <file>            : save results in BED file format of <file>
  -j  --gff <file>            : save results in GFF3 file format of <file>
  -a  --fasta <file>          : save predicted tRNA sequences in FASTA file format of <file>
  -l  --log <file>            : save log of program progress in <file>
  --detail                    : display prediction outputs in detailed view
  --brief                     : brief output format (no column headers)

  -? #                       : '#' in place of <file> chooses default name for output files
  -p  --prefix <label>        : use <label> prefix for all default output file names

  -d  --progress              : display program progress messages
  -q  --quiet                 : quiet mode (credits & run option selections suppressed)
  -y  --hitsrc                : show origin of hits (Ts=tRNAscan 1.4, Eu=EufindtRNA, 
                                  Bo=Both Ts and Eu, Inf=Infernal)

Specify Alternate Cutoffs / Data Files:

  -X  --score <score>         : set cutoff score (in bits) for reporting tRNAs (default=20)
  -g  --gencode <file>        : use alternate genetic codes specified in <file> for
                                  determining tRNA type
  -z  --pad <number>          : use <number> nucleotides padding when passing first-pass
                                  tRNA bounds predictions to CM analysis (default=8)
  --len <length>              : set max length of tRNA intron+variable region for legacy search mode
                                  (default=116bp)
Misc Options:

  -h  --help                  : print this help message
  -c  --conf <file>           : tRNAscan-SE configuration file (default: tRNAscan-SE.conf)
  -Q  --forceow               : do not prompt user before overwriting pre-existing
                                  result files  (for batch processing)

  --match <EXPR>              : search only sequences with names matching <EXPR> string
                                  (<EXPR> may contain * or ? wildcard chars)
  --search <EXPR>             : start search at sequence with name matching <EXPR> string
                                  and continue to end of input sequence file(s)
Special Advanced Options (for testing & special purposes)

  -U                          : search for tRNAs with alternate models defined in configuration file

  -t  --tscan                 : search using tRNAscan only (defaults to strict params)
  --tmode <mode>              : explicitly set tRNAscan params, where <mode>=R or S
                                  (R=relaxed, S=strict tRNAscan v1.3 params)

  -v  --verbose <file>        : save verbose tRNAscan 1.3 output to <file>
  --nomerge                   : Keep redundant tRNAscan 1.3 hits (don't filter out multiple
                                  predictions per tRNA identification)
  -e  --eufind                : search using Eukaryotic tRNA finder (EufindtRNA) only
                                  (defaults to Normal seach parameters when run alone,
                                  or to Relaxed search params when run with Cove)
  --emode <mode>              : explicitly set EufindtRNA params, where <mode>=R, N, or S
                                  (relaxed, normal, or strict)

  --iscore <score>            : manually set "intermediate" cutoff score for EufindtRNA
  -r  --fsres <file>          : save first-pass scan results from EufindtRNA, tRNAscan, or
                                  Infernal hmm in <file> in tabular results format
  --mid                       : fast scan mode - search using Infernal with mid-level strictness of hmm filter
  -F  --falsepos <file>       : save first-pass candidate tRNAs in <file> that were then
                                  found to be false positives by second-pass analysis
  --missed <file>             : save all seqs that do NOT have at least one
                                  tRNA prediction in them (aka "missed" seqs)
  --thread <number>           : number of threads used for running infernal (default is to use available threads)
```


Example

```
tRNAscan-SE \
  --thread 22 \
  -G \
  -o trnascan.res
  sample.fasta 
```