# sortmerna

主な用途は、メタトランススクリプトームデータをrRNAとそれ以外に分けること。

[github](https://github.com/biocore/sortmerna)
[github - User-manual-v4.0](https://github.com/biocore/sortmerna/wiki/User-manual-v4.0)

## インストール

```
conda install -c bioconda sortmerna
```

sortmernaオリジナルのrRNAデータベースをダウンロードする。SILVA and RFAMを合わせてrRNAじゃない配列を除いているらしい。

[new databases documented](https://github.com/biocore/sortmerna/issues/329)
[github.com/biocore/sortmerna/releases](https://github.com/biocore/sortmerna/releases)

```
cd ~/db/sortmerna/v4.3.4
wget https://github.com/biocore/sortmerna/releases/download/v4.3.4/database.tar.gz # これはv4.3.4時のデータベース（221118最新）
tar -zxvf database.tar.gz
```

解凍すると４つのファイルができる。配列数は上から73977, 31540, 309126, 225160
- smr_v4.3_default_db.fasta
- smr_v4.3_fast_db.fasta
- smr_v4.3_sensitive_db.fasta
- smr_v4.3_sensitive_db_rfam_seeds.fasta

４つのファイルの違いは、分類グループごとにどのくらいの%identityでクラスタリングしているかで、`smr_v4.3_sensitive_db.fasta`を使うのが一番センシティブ（時間はかかる）な方法だそうです（[Databases for v4.3](https://github.com/biocore/sortmerna/issues/292)）。


ちなみにsilvaのfastaをダウンロードするなら↓。rRNAじゃないものも混ざってしまっているそうなので、sortmernaオリジナルデータベースを使っておくのが無難？

```shell
cd ~/db/silva
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_LSURef_NR99_tax_silva_trunc.fasta.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta.gz
gunzip -c SILVA_138.1_LSURef_NR99_tax_silva_trunc.fasta.gz > SILVA_138.1_LSURef_NR99_tax_silva_trunc.fasta
gunzip -c SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta.gz > SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta
```



## ヘルプ

sortmerna --help

```
  -------------------------------------------------------------------------------------------------------------
  | option            type-format           description                                          default      |
  -------------------------------------------------------------------------------------------------------------

    [REQUIRED]
    --ref             PATH        Required  Reference file (FASTA) absolute or relative path.

       Use mutliple times, once per a reference file


    --reads           PATH        Required  Raw reads file (FASTA/FASTQ/FASTA.GZ/FASTQ.GZ).

       Use twice for files with paired reads.
       The file extensions are Not important. The program automatically
       recognizes the file format as flat/compressed, fasta/fastq



    [COMMON]
    --workdir         PATH        Optional  Workspace directory                         USRDIR/sortmerna/run/

       Default structure: WORKDIR/
                              idx/   (References index)
                              kvdb/  (Key-value storage for alignments)
                              out/   (processing output)
                              readb/ (pre-processed reads/index)


    --kvdb            PATH        Optional  Directory for Key-value database            WORKDIR/kvdb

       KVDB is used for storing the alignment results.


    --idx-dir         PATH        Optional  Directory for storing Reference index.      WORKDIR/idx


    --readb           PATH        Optional  Storage for pre-processed reads             WORKDIR/readb/

       Directory storing the split reads, or the random access index of compressed reads


    --fastx           BOOL        Optional  Output aligned reads into FASTA/FASTQ file
    --sam             BOOL        Optional  Output SAM alignment for aligned reads.


    --SQ              BOOL        Optional  Add SQ tags to the SAM file


    --blast           STR         Optional  output alignments in various Blast-like formats

       Sample values: '0'                    - pairwise
                      '1'                    - tabular (Blast - m 8 format)
                      '1 cigar'              - tabular + column for CIGAR
                      '1 cigar qcov'         - tabular + columns for CIGAR and query coverage
                      '1 cigar qcov qstrand' - tabular + columns for CIGAR, query coverage,
                                               and strand


    --aligned         STR/BOOL    Optional  Aligned reads file prefix [dir/][pfx]       WORKDIR/out/aligned

       Directory and file prefix for aligned output i.e. each
       output file goes into the specified directory with the given prefix.
       The appropriate extension: (fasta|fastq|blast|sam|etc) is automatically added.
       Both 'dir' and 'pfx' are optional.
       The 'dir' can be a relative or an absolute path.
       If 'dir' is not specified, the output is created in the WORKDIR/out/
       If 'pfx' is not specified, the prefix 'aligned' is used
       Examples:
       '-aligned $MYDIR/dir_1/dir_2/1' -> $MYDIR/dir_1/dir_2/1.fasta
       '-aligned dir_1/apfx'           -> $PWD/dir_1/apfx.fasta
       '-aligned dir_1/'               -> $PWD/aligned.fasta
       '-aligned apfx'                 -> $PWD/apfx.fasta
       '-aligned  (no argument)'       -> WORKDIR/out/aligned.fasta


    --other           STR/BOOL    Optional  Non-aligned reads file prefix [dir/][pfx]   WORKDIR/out/other

       Directory and file prefix for non-aligned output i.e. each
       output file goes into the specified directory with the given prefix.
       The appropriate extension: (fasta|fastq|blast|sam|etc) is automatically added.
       Must be used with 'fastx'.
       Both 'dir' and 'pfx' are optional.
       The 'dir' can be a relative or an absolute path.
       If 'dir' is not specified, the output is created in the WORKDIR/out/
       If 'pfx' is not specified, the prefix 'other' is used
       Examples:
       '-other $MYDIR/dir_1/dir_2/1' -> $MYDIR/dir_1/dir_2/1.fasta
       '-other dir_1/apfx'           -> $PWD/dir_1/apfx.fasta
       '-other dir_1/'               -> $PWD/dir_1/other.fasta
       '-other apfx'                 -> $PWD/apfx.fasta
       '-other  (no argument)'       -> aligned_out/other.fasta
                                        i.e. the same output directory
                                        as used for aligned output


    --num_alignments  INT         Optional  Positive integer (INT >=0).

       If used with '-no-best' reports first INT alignments per read reaching
       E-value threshold, which allows to lower the CPU time and memory use.
       Otherwise outputs INT best alignments.
       If INT = 0, all alignments are output


    --no-best         BOOL        Optional  Disable best alignments search                          False

       The 'best' alignment is the highest scoring alignment out of All alignments of a read,
       and the read can potentially be aligned (reaching E-value threshold) to multiple reference
       sequences.
       By default the program searches for best alignments i.e. performs an exhaustive search
       over all references. Using '-no-best' will make the program to search just
       the first N alignments, where N is set using '-num_alignments' i.e. 1 by default.


    --min_lis         INT         Optional  Search only alignments that have the LIS                2
                                            of at least N seeds long

       LIS stands for Longest Increasing Subsequence. It is computed using seeds, which
       are k-mers common to the read and the reference sequence. Sorted sequences of such seeds
       are used to filter the candidate references prior performing the Smith-Waterman alignment.


    --print_all_reads BOOL        Optional  Output null alignment strings for non-aligned reads     False
                                            to SAM and/or BLAST tabular files

    --paired          BOOL        Optional  Flags paired reads                                      False

        If a single reads file is provided, use this option to indicate
        the file contains interleaved paired reads when neither
        'paired_in' | 'paired_out' | 'out2' | 'sout' are specified.


    --paired_in       BOOL        Optional  Flags the paired-end reads as Aligned,                  False
                                            when either of them is Aligned.

        With this option both reads are output into Aligned FASTA/Q file
        Must be used with 'fastx'.
        Mutually exclusive with 'paired_out'.


    --paired_out      BOOL        Optional  Flags the paired-end reads as Non-aligned,              False
                                            when either of them is non-aligned.

        With this option both reads are output into Non-Aligned FASTA/Q file
        Must be used with 'fastx'.
        Mutually exclusive with 'paired_in'.


    --out2            BOOL        Optional  Output paired reads into separate files.                False

       Must be used with 'fastx'.
       If a single reads file is provided, this options implies interleaved paired reads
       When used with 'sout', four (4) output files for aligned reads will be generated:
       'aligned-paired-fwd, aligned-paired-rev, aligned-singleton-fwd, aligned-singleton-rev'.
       If 'other' option is also used, eight (8) output files will be generated.


    --sout            BOOL        Optional  Separate paired and singleton aligned reads.            False

       To be used with 'fastx'.
       If a single reads file is provided, this options implies interleaved paired reads
       Cannot be used with 'paired_in' | 'paired_out'


    --zip-out         STR/BOOL    Optional  Controls the output compression                       Yes/True

       By default the report files are produced in the same format as the input i.e.
       if the reads files are compressed (gz), the output is also compressed.
       The default behaviour can be overriden by using '-zip-out'.
       The possible values: Y(es), N(o), T(rue), F(alse). No value means 'True'.
       The values are Not case sensitive i.e. 'Yes, YES, yEs, Y, y' are all OK
       Examples:
       '-reads freads.gz -zip-out n' : generate flat output when the input is compressed
       '-reads freads.flat -zip-out' : compress the output when the input files are flat


    --match           INT         Optional  SW score (positive integer) for a match.                2

    --mismatch        INT         Optional  SW penalty (negative integer) for a mismatch.          -3

    --gap_open        INT         Optional  SW penalty (positive integer) for introducing a gap.    5

    --gap_ext         INT         Optional  SW penalty (positive integer) for extending a gap.      2

    -e                DOUBLE      Optional  E-value threshold.                                      1

       Defines the 'statistical significance' of a local alignment.
       Exponentially correllates with the Minimal Alignment score.
       Higher E-values (100, 1000, ...) cause More reads to Pass the alignment threshold


    -F                BOOL        Optional  Search only the forward strand.                         False

    -N                BOOL        Optional  SW penalty for ambiguous letters (N's) scored
                                            as --mismatch

    -R                BOOL        Optional  Search only the reverse-complementary strand.           False


    [OTU_PICKING]
    --id              INT         Optional  %%id similarity threshold (the alignment                0.97
                                            must still pass the E-value threshold).

    --coverage        INT         Optional  %%query coverage threshold (the alignment must          0.97
                                            still pass the E-value threshold)

    --de_novo_otu     BOOL        Optional  Output FASTA file with 'de novo' reads                  False

       Read is 'de novo' if its alignment score passes E-value threshold, but both the identity
       '-id', and the '-coverage' are below their corresponding thresholds
       i.e. ID < %%id and COV < %%cov


    --otu_map         BOOL        Optional  Output OTU map (input to QIIME's make_otu_table.py).    False
                                            Cannot be used with 'no-best because
                                            the grouping is done around the best alignment'


    [ADVANCED]
    --passes          INT,INT,INT Optional  Three intervals at which to place the seed on           L,L/2,3
                                             the read (L is the seed length)

    --edges           INT         Optional  Number (or percent if INT followed by %% sign) of       4
                                            nucleotides to add to each edge of the read
                                            prior to SW local alignment

    --num_seeds       BOOL        Optional  Number of seeds matched before searching                2
                                            for candidate LIS

    --full_search     INT         Optional  Search for all 0-error and 1-error seed                 False
                                            matches in the index rather than stopping
                                            after finding a 0-error match (<1%% gain in
                                            sensitivity with up four-fold decrease in speed)

    --pid             BOOL        Optional  Add pid to output file names.                           False

    -a                INT         Optional  DEPRECATED in favour of '-threads'. Number of           numCores
                                            processing threads to use.
                                            Automatically redirects to '-threads'

    --threads         INT         Optional  Number of Processing threads to use                     2


    [INDEXING]
    --index           INT         Optional  Build reference database index                          2

       By default when this option is not used, the program checks the reference index and
       builds it if not already existing.
       This can be changed by using '-index' as follows:
       '-index 0' - skip indexing. If the index does not exist, the program will terminate
                                and warn to build the index prior performing the alignment
       '-index 1' - only perform the indexing and terminate
       '-index 2' - the default behaviour, the same as when not using this option at all


    -L                DOUBLE      Optional  Indexing: seed length.                                  18

    -m                DOUBLE      Optional  Indexing: the amount of memory (in Mbytes) for          3072
                                            building the index.

    -v                BOOL        Optional  Produce verbose output when building the index          True

    --interval        INT         Optional  Indexing: Positive integer: index every Nth L-mer in    1
                                            the reference database e.g. '-interval 2'.

    --max_pos         INT         Optional  Indexing: maximum (integer) number of positions to      1000
                                            store for each unique L-mer.
                                            If 0 - all positions are stored.


    [HELP]
    -h                BOOL        Optional  Print help information

    --version         BOOL        Optional  Print SortMeRNA version number


    [DEVELOPER]
    --dbg_put_db      BOOL        Optional  
    --cmd             BOOL        Optional  Launch an interactive session (command prompt)          False

    --task            INT         Optional  Processing Task                                         4

       Possible values: 0 - align. Only perform alignment
                        1 - post-processing (log writing)
                        2 - generate reports
                        3 - align and post-process
                        4 - all


    --dbg-level       INT         Optional  Debug level                                             0

      Controls verbosity of the execution trace. Default value of 0 corresponds to
      the least verbose output.
      The highest value currently is 2.
```



Example

```
# when a single out of a pair of reads is aligned,
#   paired_in - put both reads into the 'aligned.fasta' file
#   paired_out - put both reads into the 'other.fasta' file

wd_rigthtend=${wd##*/} # working_dirの最小フォルダ名だけ使いたい時はこうやって取り出す

sortmerna \
 --threads 22 \
 --fastx \ # output FASTA/FASTQ file  
 --paired_in \
 --out2 \
 --zip-out \
 -v \
 --ref     ~/db/sortmerna/v4.3.4/smr_v4.3_sensitive_db.fasta \ # 複数のfastaファイルを使いたい場合は、--refをいくつも並べればok。
 --idx-dir ~/db/sortmerna/idx-dir \ # レファレンスインデックスの置き場所
 --reads   ${sample}_R1.fastq.gz \
 --reads   ${sample}_R2.fastq.gz \
 --workdir sortmerna/${sample} \
 --aligned sortmerna/${sample}/out/aligned \ # マップされた配列のfastqファイル。ファイル名はaligned_fwd.fq.gz, aligned_rev.fq.gzとなる。デフォルトだと[workdir]/out/aligned_fwd.fq.gzとかになる。
 --other   sortmerna/${sample}/out/other \ # マップされなかった配列のfastqファイル。ファイル名はother_fwd.fq.gz, other_rev.fq.gzとなる。デフォルトだとこのファイルは出ない。


# シンボリックリンクを作る場合
ln -s sortmerna/${sample}/out/other_fwd.fq.gz ${sample}_rRNAunmap_R1.fastq.gz \
ln -s sortmerna/${sample}/out/other_rev.fq.gz ${sample}_rRNAunmap_R2.fastq.gz \
```


## tRNAも探してみよう

tRNAの配列データベースを使ってsotrmeRNAをかけることで、rRNAと同じようにtRNAもフィルターしてみよう。

レファレンス用fastaは、[GtTNAdb](http://gtrnadb.ucsc.edu/)からダウンロードしたものを使ってみる。[FAQ](http://gtrnadb.ucsc.edu/faq.html)と[User Guide - How to search tRNA genes](http://gtrnadb.ucsc.edu/docs/genesearch/)に紹介されているように、tRNA Sifterから検索して検索結果のfastaファイルをダウンロードする。

2022-12-01: Domain: ALL, Minimum Score 70で検索した結果のfastaを、~/db/GtRNAdb/2022-12/に保存（270359配列）。
長さゼロの配列が入っていたようでsortmeRNAに怒られたので（seedの19baseより短い）、seqkitで除いておく。

```
$seqkit stats tRNASifter_downloaded.fa
file                      format  type  num_seqs     sum_len  min_len  avg_len  max_len
tRNASifter_downloaded.fa  FASTA   DNA    270,359  20,942,427        0     77.5      222

$seqkit seq --min-len 20 tRNASifter_downloaded.fa | seqkit stats
file  format  type  num_seqs     sum_len  min_len  avg_len  max_len
-     FASTA   DNA    269,912  20,942,427       70     77.6      222

$seqkit seq --min-len 20 tRNASifter_downloaded.fa >tRNASifter_all_minscore70.fa
```

```
sortmerna \
 --threads 22 \
 --fastx \
 --paired_in \
 --out2 \
 --zip-out \
 -v \
 --ref     ~/db/GtRNAdb/2022-12/tRNASifter_all_minscore70.fa \
 --idx-dir ~/db/sortmerna/idx-dir-trna \
 --reads   ${sample}_R1.fastq.gz \
 --reads   ${sample}_R2.fastq.gz \
 --workdir sortmerna_trna/${sample} \
 --aligned sortmerna_trna/${sample}/out/${sample}_aligned \
 --other   sortmerna_trna/${sample}/out/${sample}_other \
```

