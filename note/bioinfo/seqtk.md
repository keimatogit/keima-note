# seqtk

seqkitに似てるツール。fastqファイルでの低クオリティの箇所マスクのために入れた(220525)。

[github](https://github.com/lh3/seqtk)

```
conda install -c bioconda seqtk
```

```
Usage:   seqtk <command> <arguments>
Version: 1.3-r106

Command: seq       common transformation of FASTA/Q
         comp      get the nucleotide composition of FASTA/Q
         sample    subsample sequences
         subseq    extract subsequences from FASTA/Q
         fqchk     fastq QC (base/quality summary)
         mergepe   interleave two PE FASTA/Q files
         trimfq    trim FASTQ using the Phred algorithm

         hety      regional heterozygosity
         gc        identify high- or low-GC regions
         mutfa     point mutate FASTA at specified positions
         mergefa   merge two FASTA/Q files
         famask    apply a X-coded FASTA to a source FASTA
         dropse    drop unpaired from interleaved PE FASTA/Q
         rename    rename sequence names
         randbase  choose a random base from hets
         cutN      cut sequence at long N
         listhet   extract the position of each het
```


## seqtk seq

```
Usage:   seqtk seq [options] <in.fq>|<in.fa>

Options: -q INT    mask bases with quality lower than INT [0]
         -X INT    mask bases with quality higher than INT [255]
         -n CHAR   masked bases converted to CHAR; 0 for lowercase [0]
         -l INT    number of residues per line; 0 for 2^32-1 [0]
         -Q INT    quality shift: ASCII-INT gives base quality [33]
         -s INT    random seed (effective with -f) [11]
         -f FLOAT  sample FLOAT fraction of sequences [1]
         -M FILE   mask regions in BED or name list FILE [null]
         -L INT    drop sequences with length shorter than INT [0]
         -F CHAR   fake FASTQ quality []
         -c        mask complement region (effective with -M)
         -r        reverse complement
         -A        force FASTA output (discard quality)
         -C        drop comments at the header lines
         -N        drop sequences containing ambiguous bases
         -1        output the 2n-1 reads only
         -2        output the 2n reads only
         -V        shift quality by '(-Q) - 33'
         -U        convert all bases to uppercases
         -S        strip of white spaces in sequences
```

低クオリティの場所があるリードを削除（低クオリティの箇所をNに変換してNが含まれるリードを削除）

```
seqtk seq -q 20 -n N sample.R1.fastq >sample_maskq20.R1.fastq
seqtk seq -N sample_maskq20.R1.fastq >sample_dropq20.R1.fastq.R1.fastq
```