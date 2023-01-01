# RSEM

RSEMの中でbowtie2/STARを使ったマッピングとRSEMでのトランスクリプトのカウントを行うこともできるし、マッピング結果のBAMファイルをRSEMに渡してカウントだけしてもらうこともできる。マッピング結果のログをあまり出してくれないので、BAMを渡すことにした。

[rsem-calculate-expression.html](https://deweylab.github.io/RSEM/rsem-calculate-expression.html)

```bash
# インストール
conda install -c bioconda rsem
```
<!-- TOC depthFrom:2 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->

- [Bowtie2でマッピングして、BAMファイルをRSEMに渡す場合](#bowtie2bamrsem渡場合)
	- [1. Bowtie2でマッピング](#1-bowtie2)
	- [2. RSEMのためのインデックスを作成](#2-rsem作成)
	- [3. RSEMでカウント](#3-rsem)
- [RSEMの出力ファイル](#rsem出力)
- [（補足）出力ファイルをまとめる](#補足出力)

<!-- /TOC -->


## Bowtie2でマッピングして、BAMファイルをRSEMに渡す場合


### 1. Bowtie2でマッピング

- Bowtie2のデフォルトオプションでかけた場合RSEMが処理できなかった（何がだめだったかは控え忘れたけど、--end-to-endがないとだめなのかな？）。Trinityのalign_and_estimate_abundance.pl（bowtie2やRSEMをかけてくれるパイプライン）のヘルプの`--bowtie2_RSEM `の項で、bowtie2のオプションをどのように設定しているかが見れるので、これを参考にした。
- 渡すのは、'unsorted' BAMでないとだめなので、sortはしなくてよい。

align_and_estimate_abundance.plのbowtie2オプション設定
（align_and_estimate_abundance.plの標準出力を見ると、実際は-X 800も入ってた。-X有効なペアエンドアラインメントの最大フラグメント長で、bowtieでのデフォは500。）
```
--bowtie2_RSEM <string>         if using 'bowtie2', default: "--no-mixed --no-discordant --gbar 1000 --end-to-end -k 200 "
```



bowtie2実行とBAM作成
```shell
bowtie2 \
 --no-mixed --no-discordant --gbar 1000 --end-to-end -k 200 [-X 800] \
 -p thread \
 -x bowtie2_database \
 -1 R1.fasta \
 -2 R2.fasta \
 -S result.sam \
&& samtools view -bS result.sam > result.bam
```


### 2. RSEMのためのインデックスを作成

- fastaファイルは圧縮されてるとだめ。
- Trinityでアセンブルしたトランスクリプトへのマッピングを使う場合、`--transcript-to-gene-map`でgene名とトランスクリプト名の対応表（タブ区切り２列）を渡す。Trinity実行時に作成されたと思うけど、なければTrinityのget_Trinity_gene_to_trans_map.plでも作れるみたい。ちなみにゲノムへのマッピングを使うときは、--transcript-to-gene-mapの代わりに`--gtf`でgtfファイルを渡すそうです。gene情報を渡さなかった場合、RSEM出力でisoform_idとgene_idが同じになる。

```shell
# レファレンスがのassembled transcriptomeのとき
rsem-prepare-reference \
 -p thread \
 --transcript-to-gene-map Trinity.fasta.gene_trans_map \
 Trinity.fasta \
 rsem_index_name

# レファレンスがgenomeのとき
rsem-prepare-reference \
 -p thread \
 --gtf gene.gtf \
 genome.fasta \
 rsem_index_name
```


### 3. RSEMでカウント

```shell
rsem-calculate-expression \
 -p thread \
 --alignments --paired-end \
 bamfile.bam \
 rsem_index \
 sample_name (output_prefix)
```


## RSEMの出力ファイル

[sample_name].genes.results、[sample_name].isoforms.results、[sample_name]transcript.bam、[sample_name].stats（ディレクトリ）が作成される。genes.resultsとisoforms.resultsがそれぞれ遺伝子単位とトランスクリプト単位のカウントなどの結果（タブ区切り）。


- `transcript_id`, `gene_id`： transcript nameとそれが属するgene name。RSEMのインデックス作成時に遺伝子名情報を渡していない場合、同じ名前が入る。
- `length`: transcriptの長さ。poly-Aがついていた場合はpoly-A部分は除いた長さになるらしい。
- `effective_length`: マッピングされうる箇所の数。全長アライメントなので、[length] - [そのサンプルのmean fragment length] + 1 になるそうです（poly-Aがない場合）。mean fragment lengthがサンプルごとなので、effective_lengthもサンプルごとにちょっと違う。
- `expected_count`: 各リードがこのトランスクリプトに由来する事後確率の合計・・だそうです。似ているトランスクリプトがあって曖昧にマッピングされているリードもexpectにつかっている・・。expected_countの合計は実際のカウント合計より小さくなりがちらしい。
- `TPM`: 定義通りのTPM。長さはeffective_lengthを使用。(expected_count * 1000 / effective_length) * 1000000 / ∑(expected_count * 1000 / effective_length)
- `FPKM`: 定義通りのTPM。長さはeffective_lengthを使用。(expected_count * 1000 / effective_length) * 1000000 / ∑(expected_count)
- `IsoPct`: isoform percentage。そのisoform(transcript)が属する遺伝子の頻度に対する、そのisoform(transcript)の頻度の割合。その遺伝子のisoformがひとつのときは100。


## （補足）出力ファイルをまとめる

`rsem-generate-data-matrix`で各サンプルのisoforms/genes.resultsファイルのexpected_count列を集めて、カウントマトリックスを作れる。

```shell
# usage:
rsem-generate-data-matrix filenames > output

# example:
rsem-generate-data-matrix *.isoforms.results > all.isoforms.results
```

expected_count以外はまとめられない。自分でRとかで書いた方がいいかも。
