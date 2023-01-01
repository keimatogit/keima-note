# bcftools

インストール：
普通にやるとエラーが出るので、専用環境にしてopenssl=1.0でダウングレードする

[bcftools: error while loading shared libraries: libcrypto.so.1.0.0: cannot open shared object file: No such file or directory](https://www.biostars.org/p/9480029/)

```
conda create -n bcftools -c bioconda bcftools samtools
conda install -c bioconda openssl=1.0
```

レファレンスのfastaには.faiを作っておかないといけないみたい（`samtools faidx ref.fasta`で作れる）
出力をcompressedにするとパイプで次に渡せないので注意。
フィールドはmpileupで指定。指定なしだとGT:PLだけ出てくる。
```
samtools faidx db/Candida_albicans_NR_125332.fasta

# 1サンプルで出してみる
srun -c 8 \
samtools mpileup -Ou --max-depth 0 --output-tags AD,DP -g --reference ../db/Candida_albicans_NR_125332.fasta 03_Candida_albicans_NR_125332.sort.bam \
 | bcftools call -mv -O v -o 03.vcf

# 全サンプルで出してみる
srun -c 44 \
samtools mpileup -Ou --max-depth 0 --output-tags AD,DP -g --reference ../db/Candida_albicans_NR_125332.fasta \
 01_Candida_albicans_NR_125332.sort.bam \
 02_Candida_albicans_NR_125332.sort.bam \
 03_Candida_albicans_NR_125332.sort.bam \
 04_Candida_albicans_NR_125332.sort.bam \
 05_Candida_albicans_NR_125332.sort.bam \
 06_Candida_albicans_NR_125332.sort.bam \
 07_Candida_albicans_NR_125332.sort.bam \
 08_Candida_albicans_NR_125332.sort.bam \
 09_Candida_albicans_NR_125332.sort.bam \
 NTC_Candida_albicans_NR_125332.sort.bam \
 | bcftools call -mv -O v -o all_Candida_albicans_NR_125332.vcf

```

```
>samtools mpileup -h
mpileup: option requires an argument -- 'h'

Usage: samtools mpileup [options] in1.bam [in2.bam [...]]

Input options:
  -6, --illumina1.3+      quality is in the Illumina-1.3+ encoding
  -A, --count-orphans     do not discard anomalous read pairs
  -b, --bam-list FILE     list of input BAM filenames, one per line
  -B, --no-BAQ            disable BAQ (per-Base Alignment Quality)
  -C, --adjust-MQ INT     adjust mapping quality; recommended:50, disable:0 [0]
  -d, --max-depth INT     max per-file depth; avoids excessive memory usage [250]
  -E, --redo-BAQ          recalculate BAQ on the fly, ignore existing BQs
  -f, --fasta-ref FILE    faidx indexed reference sequence file
  -G, --exclude-RG FILE   exclude read groups listed in FILE
  -l, --positions FILE    skip unlisted positions (chr pos) or regions (BED)
  -q, --min-MQ INT        skip alignments with mapQ smaller than INT [0]
  -Q, --min-BQ INT        skip bases with baseQ/BAQ smaller than INT [13]
  -r, --region REG        region in which pileup is generated
  -R, --ignore-RG         ignore RG tags (one BAM = one sample)
  --rf, --incl-flags STR|INT  required flags: skip reads with mask bits unset []
  --ff, --excl-flags STR|INT  filter flags: skip reads with mask bits set
                                            [UNMAP,SECONDARY,QCFAIL,DUP]
  -x, --ignore-overlaps   disable read-pair overlap detection

Output options:
  -o, --output FILE       write output to FILE [standard output]
  -g, --BCF               generate genotype likelihoods in BCF format
  -v, --VCF               generate genotype likelihoods in VCF format

Output options for mpileup format (without -g/-v):
  -O, --output-BP         output base positions on reads
  -s, --output-MQ         output mapping quality

Output options for genotype likelihoods (when -g/-v is used):
  -t, --output-tags LIST  optional tags to output:
               DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR []
  -u, --uncompressed      generate uncompressed VCF/BCF output

SNP/INDEL genotype likelihoods options (effective with -g/-v):
  -e, --ext-prob INT      Phred-scaled gap extension seq error probability [20]
  -F, --gap-frac FLOAT    minimum fraction of gapped reads [0.002]
  -h, --tandem-qual INT   coefficient for homopolymer errors [100]
  -I, --skip-indels       do not perform indel calling
  -L, --max-idepth INT    maximum per-file depth for INDEL calling [250]
  -m, --min-ireads INT    minimum number gapped reads for indel candidates [1]
  -o, --open-prob INT     Phred-scaled gap open seq error probability [40]
  -p, --per-sample-mF     apply -m and -F per-sample for increased sensitivity
  -P, --platforms STR     comma separated list of platforms for indels [all]
      --input-fmt-option OPT[=VAL]
               Specify a single input file format option in the form
               of OPTION or OPTION=VALUE
      --reference FILE
               Reference sequence FASTA FILE [null]

Notes: Assuming diploid individuals.
```


```
>bcftools call -h

About:   SNP/indel variant calling from VCF/BCF. To be used in conjunction with samtools mpileup.
         This command replaces the former "bcftools view" caller. Some of the original
         functionality has been temporarily lost in the process of transition to htslib,
         but will be added back on popular demand. The original calling model can be
         invoked with the -c option.
Usage:   bcftools call [options] <in.vcf.gz>

File format options:
       --no-version                do not append version and command line to the header
   -o, --output <file>             write output to a file [standard output]
   -O, --output-type <b|u|z|v>     output type: 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
       --ploidy <assembly>[?]      predefined ploidy, 'list' to print available settings, append '?' for details
       --ploidy-file <file>        space/tab-delimited list of CHROM,FROM,TO,SEX,PLOIDY
   -r, --regions <region>          restrict to comma-separated list of regions
   -R, --regions-file <file>       restrict to regions listed in a file
   -s, --samples <list>            list of samples to include [all samples]
   -S, --samples-file <file>       PED file or a file with an optional column with sex (see man page for details) [all samples]
   -t, --targets <region>          similar to -r but streams rather than index-jumps
   -T, --targets-file <file>       similar to -R but streams rather than index-jumps
       --threads <int>             number of extra output compression threads [0]

Input/output options:
   -A, --keep-alts                 keep all possible alternate alleles at variant sites
   -f, --format-fields <list>      output format fields: GQ,GP (lowercase allowed) []
   -g, --gvcf <int>,[...]          group non-variant sites into gVCF blocks by minimum per-sample DP
   -i, --insert-missed             output also sites missed by mpileup but present in -T
   -M, --keep-masked-ref           keep sites with masked reference allele (REF=N)
   -V, --skip-variants <type>      skip indels/snps
   -v, --variants-only             output variant sites only

Consensus/variant calling options:
   -c, --consensus-caller          the original calling method (conflicts with -m)
   -C, --constrain <str>           one of: alleles, trio (see manual)
   -m, --multiallelic-caller       alternative model for multiallelic and rare-variant calling (conflicts with -c)
   -n, --novel-rate <float>,[...]  likelihood of novel mutation for constrained trio calling, see man page for details [1e-8,1e-9,1e-9]
   -p, --pval-threshold <float>    variant if P(ref|D)<FLOAT with -c [0.5]
   -P, --prior <float>             mutation rate (use bigger for greater sensitivity) [1.1e-3]

```