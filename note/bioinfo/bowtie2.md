# Bowtie2

老舗マッピングツール

- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
- [bioinformatics/Bowtie2](https://bi.biopapyrus.jp/rnaseq/mapping/bowtie2/)

インストール
```
conda install -c bioconda bowtie2
```


## インデックスの準備

fastaから作成（.gz可）
```sh
bowtie2-build --threads 22 fasta(.gz) index_name 
```

または、[Bowtie2ホームページ](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)のIndex欄からダウンロード -> unzip。



## マッピング

- `-x`：bowtie2のインデックス
- `-1`：ペアエンドのR1ファイル（.gz可）
- `-2`：ペアエンドのR2ファイル（.gz可）
- `-p`：スレッド数
- `-S`：マッピング結果のsamファイル（指定がない場合は標準出力）

```sh
bowtie2 \
 -x bowtie2_index \
 -1 sample_R1.fastq.gz \
 -2 sample_R2.fastq.gz \
 -p 22 \
 -b sample.bam \
 -S sample.sam
```


## RSEM用マッピング

Trinityのalign_and_estimate_abundance.plでbowtie2 -> RSEMができる。RSEMでもbowtie2 -> RSEMができる。bowtie2を個別でかけて結果のbamファイルをRSEMにかける場合は、`align_and_estimate_abundance.pl`のbowtie2のデフルトォオプション（以下）を参照（実際はこれに加えて-X 800が入っている）。RSEMには"unsorted" bamが必要。なので作成しておく。

```sh
--bowtie2_RSEM <string>         if using 'bowtie2', default: "--no-mixed --no-discordant --gbar 1000 --end-to-end -k 200 "
```

マッピング
```sh
bowtie2 \
 --no-mixed --no-discordant --gbar 1000 --end-to-end -k 200 [-X 800]\
 -x bowtie2_index \
 -1 R1.fastq.gz \
 -2 R2.fastq.gz \
 -p 22 \
 -S samfile.sam \
&& samtools view -bS samfile.sam > bamfile.bam
```

