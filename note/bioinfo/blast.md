# blast

[Heavy Watal - BLAST](https://heavywatal.github.io/bio/blast.html)

```
conda install -c bioconda blast
```



makeblastdb
元になるfastaは圧縮されているとだめなので、gunzip -cで標準入力から入れる（このとき -titleと-outが必須）。
```
gunzip -c KEGG.ftp_210909/fasta/eukaryotes.pep.gz | \
  makeblastdb -in - -title eukaryotes -out blast/kegg/eukaryotes -dbtype prot
```

例）nrのblastdbを作る
```
# get fasta
srun wget https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz

# make blast db
srun -c 44 gunzip -c ~/db/ncbi_fasta/nr.gz | \
  makeblastdb -in - -title nr -out ~/db/blast/nr -dbtype prot -parse_seqids
> BLAST Database creation error: Defline lacks a proper ID around line 414241776
```

nt

```
sbatch -c 22 -p c6420 -J wget --wrap="wget https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz"
sbatch -J makeblastdb -p c6420 -c 44 --wrap="gunzip -c ~/db/ncbi/fasta/nt.gz | makeblastdb -in - -title nt -out ~/db/blast/nt -dbtype nucl -parse_seqids"

# エラーになってしまった
BLAST Database creation error: Error: Duplicate seq_ids are found: 
LCL|6O9K_A
```

Example
-queryは圧縮ファイルはだめかも。fastqもだめかも。fastaのみ？
```
blastp \
 -num_threads 22 \
 -outfmt 6 \
 -query query.fasta \
 -out   query.nr.blast \
 -db    ~/db/blast/nr  \
 -max_target_seqs 100 \
 -evalue 1e-5
```


outfmt 6 (タブ区切り)のデフォ
1.  qseqid      query or source (e.g., gene) sequence id
2.  sseqid      subject  or target (e.g., reference genome) sequence id
3.  pident      percentage of identical matches
4.  length      alignment length (sequence overlap)
5.  mismatch    number of mismatches
6.  gapopen     number of gap openings
7.  qstart      start of alignment in query
8.  qend        end of alignment in query
9.  sstart      start of alignment in subject
10.  send        end of alignment in subject
11.  evalue      expect value
12.  bitscore    bit score