# abricate

[abricate](https://github.com/tseemann/abricate)

```
conda create -n abricate
conda install -c conda-forge -c bioconda -c defaults abricate
```

手持ちのデータベース一覧
```
abricate --list
```




```shell
source .bashrc.imetconda
conda activate abricate
cd ~/211013_abricate/0_fastq

# Symbolic link to fastq
ln -s /imetgpfs/projects/cw/IMET/211008-IMET79-c049/02_iMetDBnt-bwa/*hg19-bwa.unmapped.fastq.gz .

# Assemble
for i in IMET0923 IMET0924
do sbatch \
  -c 96 \
  -J megahit_${i} \
  -o slurm-megahit-${i}-%j.out \
  --export=\
c=96,\
R1=./0_fastq/${i}.R1.hg19-bwa.unmapped.fastq.gz,\
R2=./0_fastq/${i}.R2.hg19-bwa.unmapped.fastq.gz,\
outdir=./1_megahit/${i} \
  ~/sh/megahit.sh
done
# symbolic link
for i in IMET0923 IMET0924
do ln -s \
  ${i}/final.contigs.fa.gz contigs_${i}.fasta.gz
done

# Symbolic link to contigs of TKDA1138(コロナhuman gutショットガン)
# moderate
ln -s /imetgpfs/projects/cw/TKDA/Maeda/210203-TKDA1138-v192/04_gene_content/1_megahit/contigs_20[456].fasta.gz .
# control
ln -s /imetgpfs/projects/cw/TKDA/Maeda/210203-TKDA1138-v192/04_gene_content/1_megahit/contigs_30[123].fasta.gz .
```
実行テスト
```shell
srun -c 96 \
abricate \
 1_megahit/contigs_20[456].fasta.gz \
 --threads 96 > 2_abricate/204_5_6
```

~/sh/abricate.sh
```bash
#!/bin/bash
#SBATCH -c 22

source /imetgpfs/miniconda3/bin/activate abricate

abricate \
 ${in} \
 --db ${db} \
 --threads 22 > ${out}
```


```shell
# TKDAサンプル
for i in 204 205 206 301 302 303
do sbatch \
  -J abricate_resfinder_${i} \
  -o slurm-abricate-resfinder-${i}-%j.out \
  --export=\
in=./1_megahit/contigs_${i}.fasta.gz,\
db=resfinder,\
out=./2_abricate/abricate_resfinder_${i}.tsv\
  ~/sh/abricate.sh
done

# IMETサンプル
for i in IMET0923 IMET0924
do sbatch \
  -J abricate_resfinder_${i} \
  -o slurm-abricate-resfinder-${i}-%j.out \
  --export=\
in=./1_megahit/contigs_${i}.fasta.gz,\
db=ncbi,\
out=./2_abricate/abricate_ncbi_${i}.tsv\
  ~/sh/abricate.sh
done

# まとめる
srun abricate --summary abricate_ncbi_* > summary_ncbi.tsv
srun abricate --summary abricate_resfinder_* > summary_resfinder.tsv

# まとめてかける
sbatch \
  -J abricate_vfdb \
  -o slurm-abricate-vfdb-%j.out \
  --export=\
in=./1_megahit/contigs_*.fasta.gz,\
db=vfdb,\
out=./2_abricate/abricate_vfdb_all.tsv\
  ~/sh/abricate.sh
srun abricate --summary abricate_vfdb_all.tsv > summary_vfdb.tsv
```

元岡さんに報告して、/imetgpfs/projects/cw/TKDA/Maeda/210203-TKDA1138-v192の健常サンプルは全部かけることになった。ついでに他のもかけてしまおうかな。
```shell
cd /imetgpfs/projects/cw/TKDA/Maeda/210203-TKDA1138-v192/04_gene_content

for i in 0 1 2 3
do sbatch \
  -J abricate_card_${i}xx \
  -o slurm-abricate_card_${i}xx-%j.out \
  --export=\
in=./1_megahit/contigs_${i}*.fasta.gz,\
db=card,\
out=./2_abricate/card_${i}xx.tsv\
  ~/sh/abricate.sh
done

cat ncbi* > ncbi_all.tsv
srun abricate --summary ncbi_all.tsv > summary_ncbi_all.tsv
```
