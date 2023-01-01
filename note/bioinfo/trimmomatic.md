# MEGAHIT

- [Trimmomatic Manual: V0.32 - USADELLAB.org(PDF)](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf)
- [bioinformatics/Trimmomatic](https://bi.biopapyrus.jp/rnaseq/qc/trimmomatic.html)

インストール
```
conda install -c bioconda trimmomatic
```

## ヘルプ

```
Usage: 
       PE [-version] [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-summary <statsSummaryFile>] [-quiet] [-validatePairs] [-basein <inputBase> | <inputFile1> <inputFile2>] [-baseout <outputBase> | <outputFile1P> <outputFile1U> <outputFile2P> <outputFile2U>] <trimmer1>...
   or: 
       SE [-version] [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-summary <statsSummaryFile>] [-quiet] <inputFile> <outputFile> <trimmer1>...
   or: 
       -version
```


## Example

Paired end
```
trimmomatic PE \
 -phred33 \
 -threads 22 \
 ${sample}_R1.fastq.gz       ${sample}_R2.fastq.gz \
 ${sample}_R1_trimmed.fastq  ${sample}_R1_trimmed.unpaired.fastq \
 ${sample}_R2_trimmed.fastq  ${sample}_R2_trimmed.unpaired.fastq \
 ILLUMINACLIP:adapters.fasta:2:30:10 \
 LEADING:20 \ # リードの先頭からクオリティスコアが 20 未満の塩基を切り捨てる。
 TRAILING:20 \ # リードの後尾からクオリティスコアが 20 未満の塩基を切り捨てる。
 SLIDINGWINDOW:4:15 \ # 最初の 4 はウィンドウサイズ、最後の 15 は平均クオリティ。ウィンドウ内部の平均クオリティが指定値よりも低ければ以降のリードを切り捨てる。
 MINLEN:100


# 必要に応じて
gzip ${sample}_R1_trimmed.fastq ${sample}_R2_trimmed.fastq \
rm ${sample}_R1_trimmed.unpaired.fastq ${sample}_R2_trimmed.unpaired.fastq \
&& fastqc -t 22 -o fastqc ${sample}_R1_trimmed.fastq ${sample}_R2_trimmed.fastq
```

`-phred33`はつけなくても変わらなかった（デフォでphred33）けど、心配ならつけとくといい。