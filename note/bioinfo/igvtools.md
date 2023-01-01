# igvtools

[Running igvtools from the Command Line](https://software.broadinstitute.org/software/igv/igvtools_commandline)
[igvtoolsはすぐれものーリードをカウントして各塩基の塩基組成のファイルを作る](https://ncrna.jp/blog/item/65-igv-toolkit)

```
conda create -n igvtools
conda activate igvtools

conda install -c bioconda igvtools igv
conda install -c conda-forge openjdk
```

1塩基ごとに頻度をカウント。outputのファイル名は.wig(テキスト)か.tdf（バイナリ？）かstdinのみ。
```
do igvtools count \
 -w 1 \
 --bases  \
 sample.bam \
 sample_count.wig \
 reference.fasta
```