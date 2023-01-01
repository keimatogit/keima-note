# Diamond blast

blastpとblastxをすごく早く行う。メタゲでよく使われている。

[GitHub Wiki](https://github.com/bbuchfink/diamond/wiki)
[macでインフォマティクス](https://kazumaxneo.hatenablog.com/entry/2017/08/24/150512)
[Buchfink et al. (2015). Fast and sensitive protein alignment using DIAMOND. Nature methods, 12(1), 59-60.](https://www.nature.com/articles/nmeth.3176)

```shell
# インストール
conda install -c bioconda diamond
```

##

インデックス作成。配列ファイルはgzip可、複数可。--inがない場合は標準入力を使用。
```
srun diamond makedb --in amino_acids.fasta -d database_name [--taxonmap ] [--taxonnodes
```

複数ファイルからDB
```
 zcat KEGG.ftp_210909/fasta/prokaryotes.pep.gz \
      KEGG.ftp_210909/fasta/eukaryotes.pep.gz \
      GenomeNet.ftp/mgenes/meta.pep.gz |
 diamond makedb \
  --threads 96 \
  -d diamond/kegg/prokaryotes_eukaryotes_mgenesMeta
```

nr

```
wget https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
diamond makedb --in nr.gz -d diamond/nr
```

CAZy

```
wget https://bcb.unl.edu/dbCAN2/download/CAZyDB.08062022.fa
diamond makedb --in ~/db/dbCAN/CAZyDB.08062022.fa -d ~/db/diamond/CAZyDB.08062022
```

## Example

```
diamond blastp \
  --threads 22 \
  --outfmt 6 \
  --min-score 60 \
  --query-cover 50 \
  -q input_aa.fasta \
  -d ~/db/diamond/CAZyDB.08062022 \
  -o input_aa.daa
  
```

## help

```
diamond help
diamond v0.9.14.115 | by Benjamin Buchfink <buchfink@gmail.com>
Licensed under the GNU AGPL <https://www.gnu.org/licenses/agpl.txt>
Check http://github.com/bbuchfink/diamond for updates.

Syntax: diamond COMMAND [OPTIONS]

Commands:
makedb	Build DIAMOND database from a FASTA file
blastp	Align amino acid query sequences against a protein reference database
blastx	Align DNA query sequences against a protein reference database
view	View DIAMOND alignment archive (DAA) formatted file
help	Produce help message
version	Display version information
getseq	Retrieve sequences from a DIAMOND database file
dbinfo	Print information about a DIAMOND database file

General options:
--threads (-p)         number of CPU threads
--db (-d)              database file
--out (-o)             output file
--outfmt (-f)          output format
	0   = BLAST pairwise
	5   = BLAST XML
	6   = BLAST tabular
	100 = DIAMOND alignment archive (DAA)
	101 = SAM

	Value 6 may be followed by a space-separated list of these keywords:

	qseqid means Query Seq - id
	qlen means Query sequence length
	sseqid means Subject Seq - id
	sallseqid means All subject Seq - id(s), separated by a ';'
	slen means Subject sequence length
	qstart means Start of alignment in query
	qend means End of alignment in query
	sstart means Start of alignment in subject
	send means End of alignment in subject
	qseq means Aligned part of query sequence
	sseq means Aligned part of subject sequence
	evalue means Expect value
	bitscore means Bit score
	score means Raw score
	length means Alignment length
	pident means Percentage of identical matches
	nident means Number of identical matches
	mismatch means Number of mismatches
	positive means Number of positive - scoring matches
	gapopen means Number of gap openings
	gaps means Total number of gaps
	ppos means Percentage of positive - scoring matches
	qframe means Query frame
	btop means Blast traceback operations(BTOP)
	staxids means unique Subject Taxonomy ID(s), separated by a ';' (in numerical order)
	stitle means Subject Title
	salltitles means All Subject Title(s), separated by a '<>'
	qcovhsp means Query Coverage Per HSP
	qtitle means Query title

	Default: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
--verbose (-v)         verbose console output
--log                  enable debug log
--quiet                disable console output

Makedb options:
--in                   input reference file in FASTA format

Aligner options:
--query (-q)           input query file
--strand               query strands to search (both/minus/plus)
--un                   file for unaligned queries
--unal                 report unaligned queries (0=no, 1=yes)
--max-target-seqs (-k) maximum number of target sequences to report alignments for
--top                  report alignments within this percentage range of top alignment score (overrides --max-target-seqs)
--overlap-culling      delete hits only if a higher scoring hit envelops the given percentage of the query range
--compress             compression for output files (0=none, 1=gzip)
--evalue (-e)          maximum e-value to report alignments (default=0.001)
--min-score            minimum bit score to report alignments (overrides e-value setting)
--id                   minimum identity% to report an alignment
--query-cover          minimum query cover% to report an alignment
--subject-cover        minimum subject cover% to report an alignment
--sensitive            enable sensitive mode (default: fast)
--more-sensitive       enable more sensitive mode (default: fast)
--block-size (-b)      sequence block size in billions of letters (default=2.0)
--index-chunks (-c)    number of chunks for index processing
--tmpdir (-t)          directory for temporary files
--gapopen              gap open penalty
--gapextend            gap extension penalty
--frameshift (-F)      frame shift penalty (default=disabled)
--matrix               score matrix for protein alignment (default=BLOSUM62)
--custom-matrix        file containing custom scoring matrix
--lambda               lambda parameter for custom matrix
--K                    K parameter for custom matrix
--comp-based-stats     enable composition based statistics (0/1=default)
--masking              enable masking of low complexity regions (0/1=default)
--query-gencode        genetic code to use to translate query (see user manual)
--salltitles           include full subject titles in DAA file
--sallseqid            include all subject ids in DAA file
--no-self-hits         suppress reporting of identical self hits
--taxonmap             protein accession to taxid mapping file
--taxonnodes           taxonomy nodes.dmp from NCBI

Advanced options:
--algo                 Seed search algorithm (0=double-indexed/1=query-indexed)
--bin                  number of query bins for seed search
--min-orf (-l)         ignore translated sequences without an open reading frame of at least this length
--freq-sd              number of standard deviations for ignoring frequent seeds
--id2                  minimum number of identities for stage 1 hit
--window (-w)          window size for local hit search
--xdrop (-x)           xdrop for ungapped alignment
--ungapped-score       minimum alignment score to continue local extension
--hit-band             band for hit verification
--hit-score            minimum score to keep a tentative alignment
--gapped-xdrop (-X)    xdrop for gapped alignment in bits
--band                 band for dynamic programming computation
--shapes (-s)          number of seed shapes (0 = all available)
--shape-mask           seed shapes
--index-mode           index mode (0=4x12, 1=16x9)
--rank-ratio           include subjects within this ratio of last hit (stage 1)
--rank-ratio2          include subjects within this ratio of last hit (stage 2)
--max-hsps             maximum number of HSPs per subject sequence to save for each query
--dbsize               effective database size (in letters)
--no-auto-append       disable auto appending of DAA and DMND file extensions
--xml-blord-format     Use gnl|BL_ORD_ID| style format in XML output

View options:
--daa (-a)             DIAMOND alignment archive (DAA) file
--forwardonly          only show alignments of forward strand

Getseq options:
--seq                  Sequence numbers to display.
```