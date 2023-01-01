# itsx

[ITSx - Microbiology.se](https://microbiology.se/software/itsx/)
[User's guide:(PDF)](https://microbiology.se/publ/itsx_users_guide.pdf)
[Bengtsson‐Palme, J. et al](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12073)

```
conda install -c bioconda itsx
```

```
> ITSx --help

ITSx -- Identifies ITS sequences and extracts the ITS region
by Johan Bengtsson-Palme et al., University of Gothenburg
Version: 1.1.3
-----------------------------------------------------------------
Usage: ITSx -i <input file> -o <output file>
Options:
-i {file} : DNA FASTA input file to investigate
-o {file} : Base for the names of output file(s)
-p {directory} : A path to a directory of HMM-profile collections representing ITS conserved regions, default is in the same directory as ITSx itself
--stdin {T or F} : Use input from standard input instead of an input file, off (F) by default
--date {T or F} : Adds a date and time stamp to the output directory, off (F) by default
--reset {T or F} : Re-creates the HMM-database before ITSx is run, off (F) by default

Sequence selection options:
-t {character code} : Profile set to use for the search, see the User's Guide (comma-separated), default is all
-E {value} : Domain E-value cutoff for a sequence to be included in the output, default = 1e-5
-S {value} : Domain score cutoff for a sequence to be included in the output, default = 0
-N {value} : The minimal number of domains that must match a sequence before it is included, default = 2
--selection_priority {score, sum, domains, eval} : Selects what will be of highest priority when determining the origin of the sequence, default is score
--search_eval {value} : The E-value cutoff used in the HMMER search, high numbers may slow down the process, cannot be used with the --search_score option, default is to use score cutoff, not E-value
--search_score {value} : The score cutoff used in the HMMER search, low numbers may slow down the process, cannot be used with the --search_eval option, default = 0
--allow_single_domain {e-value,score or F} : Allow inclusion of sequences that only find a single domain, given that they meet the given E-value and score thresholds, on with parameters 1e-9,0 by default
--allow_reorder {T or F} : Allows profiles to be in the wrong order on extracted sequences, off (F) by default
--complement {T or F} : Checks both DNA strands against the database, creating reverse complements, on (T) by default
--cpu {value} : the number of CPU threads to use, default is 1
--multi_thread {T or F} : Multi-thread the HMMER-search, on (T) if number of CPUs (--cpu option > 1), else off (F) by default
--heuristics {T or F} : Selects whether to use HMMER's heuristic filtering, off (F) by default
--nhmmer {T or F} : Selects whether to use nhmmer instead of hmmsearch for HMMER searches, off (F) by default

Output options:
--summary {T or F} : Summary of results output, on (T) by default
--graphical {T or F} : 'Graphical' output, on (T) by default
--fasta {T or F} : FASTA-format output of extracted ITS sequences, on (T) by default
--preserve {T or F} : Preserve sequence headers in input file instead of printing out ITSx headers, off (F) by default
--save_regions {SSU,ITS1,5.8S,ITS2,LSU,all,none} : A comma separated list of regions to output separate FASTA files for, 'ITS1,ITS2' by default
--anchor {integer or HMM} : Saves an additional number of bases before and after each extracted region. If set to 'HMM' all bases matching the corresponding HMM will be output, default = 0
--require_anchor {integer or HMM} : Requires the complete anchor to found in order to be included in the output sequences (see --anchor above). Cannot be used together with the --anchor option, default = 0
--only_full {T or F} : If true, output is limited to full-length regions, off (F) by default
--partial {integer} : Saves additional FASTA-files for full and partial ITS sequences longer than the specified cutoff, default = 0 (off)
--concat {T or F} : Saves a FASTA-file with concatenated ITS sequences (with 5.8S removed), off (F) by default
--minlen {integer} : Minimum length the ITS regions must be to be outputted in the concatenated file (see above), default = 0
--positions {T or F} : Table format output containing the positions ITS sequences were found in, on (T) by default
--table {T or F} : Table format output of sequences containing probable ITS sequences, off (F) by default
--not_found {T or F} : Saves a list of non-found entries, on (T) by default
--detailed_results {T or F} : Saves a tab-separated list of all results, off (F) by default
--truncate {T or F} : Truncates the FASTA output to only contain the actual ITS sequences found, on (T) by default
--silent {T or F} : Supresses printing progress info to stderr, off (F) by default
--graph_scale {value} : Sets the scale of the graph output, if value is zero, a percentage view is shown, default = 0
--save_raw {T or F} : Saves all raw data for searches etc. instead of removing it on finish, off (F) by default
--temp {directory} : Custom directory to put the temporary files in

-h : displays this help message
--help : displays this help message
--bugs : displays the bug fixes and known bugs in this version of ITSx
--license : displays licensing information
-----------------------------------------------------------------
```

```
sbatch -o slurm_itsx_%j.out -c 22 -p c6420 --wrap="ITSx --cpu 22 --multi_thread T --graphical F --save_regions ITS1 --positions F -i NKTAall_R1_head.fasta -o NKTAall_R1_head_ITSx"
```

fastqは扱えないみたい・・・-> ITSxpressというやつは扱えるらしい。
結構時間がかかる。fungiのアンプリコン10000リードで3分ちょっと、992228リードだと8時間半くらいかかりそう（見込み）。ITSxpressの方がだいぶ早くはある。

[darcyj/fastq-from-ITSx](https://github.com/darcyj/fastq-from-ITSx)：fastqを扱えるようなものを作った人もいるらしいけど、この人はITSxpressをお勧めしてるみたい。