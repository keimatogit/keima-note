# QIIME2

- [QIIM2docs](https://docs.qiime2.org/2022.8/)
- [QIIM2docs - “Moving Pictures” tutorial](https://docs.qiime2.org/2022.8/tutorials/moving-pictures)
- [github - otagoedna/getting_started_with_qiime2](https://github.com/otagoedna/getting_started_with_qiime2)
- [Qiime2 を用いた 16S rRNA 菌叢解析](https://qiita.com/keisuke-ota/items/6399b2f2f7459cd9e418)

短い16Sを使用する場合、unifrac距離用の系統樹はSEPPというのを使ったq2-fragment-insertionを使うといいらしい（既存treeの中に配列を配置する）。既存treeの配列と75%以上の総同性がない配列は捨てられる。
- [Ranking the biases: The choice of OTUs vs. ASVs in 16S rRNA amplicon data analysis has stronger effects on diversity measures than rarefaction and OTU identity threshold](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0264443#sec002)
- [Plugin website](https://github.com/qiime2/q2-fragment-insertion)
- [q2-fragment-insertion](https://library.qiime2.org/plugins/q2-fragment-insertion/16/)


QIIME2とphyloseqのUnifrac距離は随分違うらしい
- [very different weighted unifrac values for qiime2 versus phyloseq](https://github.com/joey711/phyloseq/issues/956)

qiime2Rの参考サイト
- [Qiime2R Github page](https://github.com/jbisanz/qiime2R)
- [Qiime2R Tutorial in Qiime2 Forum](https://forum.qiime2.org/t/tutorial-integrating-qiime2-and-r-for-data-visualization-and-analysis-using-qiime2r/4121)

- ファイルは.qza形式（zip的なもの？）で取り扱う。
- .qzvは[QIIME2 view](https://view.qiime2.org/)で見れる形式らしい。
- Rではqzaファイルを扱うqiime2Rというライブラリがある。


## installing (Miniconda)

[Natively installing QIIME 2 - Miniconda](https://docs.qiime2.org/2022.2/install/native/#miniconda)

```
wget https://data.qiime2.org/distro/core/qiime2-2022.2-py38-linux-conda.yml
conda env create -n qiime2-2022.2 --file qiime2-2022.2-py38-linux-conda.yml
# OPTIONAL CLEANUP
rm qiime2-2022.2-py38-linux-conda.yml
```

2022.8バージョン

```
wget https://data.qiime2.org/distro/core/qiime2-2022.8-py38-linux-conda.yml
conda env create -n qiime2-2022.8 --file qiime2-2022.8-py38-linux-conda.yml
# OPTIONAL CLEANUP
rm qiime2-2022.8-py38-linux-conda.yml
```

アノテーション用のデータベースは、[QIIME2 docs - Data resources](https://docs.qiime2.org/2022.8/data-resources)からダウンロード。

```
cd ~/db/qiime2

# Naive Bayes用の学習済み分類器
# （Taxonomy classifiers for use with q2-feature-classifier - Weighted Taxonomic Classifiers）
wget https://data.qiime2.org/2022.8/common/silva-138-99-nb-classifier.qza

# Blast・Vsearch用の配列と分類群
# （Marker gene reference databases - Silva (16S/18S rRNA)）
wget https://data.qiime2.org/2022.8/common/silva-138-99-seqs.qza
wget https://data.qiime2.org/2022.8/common/silva-138-99-tax.qza
```

Rパッケージqiime2Rはgithubからダウンロード（condaにはまだないみたい）

```R
devtools::install_github("jbisanz/qiime2R")
```

## Example: テーブル作成まで

[github - otagoedna/getting_started_with_qiime2 - Example Workflow](https://github.com/otagoedna/getting_started_with_qiime2/blob/master/first_workflow.md)
[Overview of QIIME 2 Plugin Workflows](https://docs.qiime2.org/2022.8/tutorials/overview/#)
[Qiime2 を用いた 16S rRNA 菌叢解析](https://qiita.com/keisuke-ota/items/6399b2f2f7459cd9e418)


### 先にプライマー配列を除いておいた方がいいかも

```
cutadapt \
 -g AGAGTTTGATYMTGGCTCAG \
 -G GCTGCCTCCCGTAGGAGT \
 --quality-cutoff 20 \
 --minimum-length 100 \
 -o ${sample}_cutadapt_R1.fastq.gz \
 -p ${sample}_cutadapt_R2.fastq.gz \
 ${sample}_R1.fastq.gz \
 ${sample}_R2.fastq.gz
```

インポート後にqiime cutadapt trim-pairedでcutadaptをかけれるみたい。


### 0) サンプルメタデータの作成

初めにサンプル情報のtsvを作成しておこう。1列目はsampleid（import時のマニフェストファイルの一番左と同じ）で、続く列にサンプル情報を格納しよう。

ここではメタデータファイルを${prefix}_metadata.tsvとする（後々使う）。

### 1) import - FastqデータをインポートしてQIIME2形式にする

[Importing data](https://docs.qiime2.org/2022.8/tutorials/importing/)

Demultiplesされたファイルをインポートする場合は、まず「マニフェストファイル」を作成する。マニフェストファイルはタブ区切りで、１列目：サンプルid、２列目：R1のフルパス、3列目：R2のフルパス（$HOMEや$PWDも使える）。

```
sample-id     forward-absolute-filepath       reverse-absolute-filepath
sample-1      $PWD/some/filepath/sample1_R1.fastq.gz  $PWD/some/filepath/sample1_R2.fastq.gz
sample-2      $PWD/some/filepath/sample2_R1.fastq.gz  $PWD/some/filepath/sample2_R2.fastq.gz
sample-3      $PWD/some/filepath/sample3_R1.fastq.gz  $PWD/some/filepath/sample3_R2.fastq.gz
```

シングルリードの場合は２列。

```
sample-id     absolute-filepath
sample-1      $PWD/some/filepath/sample1_R1.fastq
sample-2      $PWD/some/filepath/sample2_R1.fastq
```

マニフェストファイルを使ってimport。パスが間違っていても何も警告してくれないみたいなので注意。

```
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.txt \
  --input-format PairedEndFastqManifestPhred33V2 \
  --output-path ~/mykinso/qiime2/${prefix}.qza
```

qiime tools importのヘルプ。--typeや--input-formatで設定出来る値は、--show-importable-typesや--show-importable-formatsでチェック。

```
Usage: qiime tools import [OPTIONS]

  Import data to create a new QIIME 2 Artifact. See https://docs.qiime2.org/
  for usage examples and details on the file types and associated semantic
  types that can be imported.

Options:
  --type TEXT             The semantic type of the artifact that will be
                          created upon importing. Use --show-importable-types
                          to see what importable semantic types are available
                          in the current deployment.                [required]
  --input-path PATH       Path to file or directory that should be imported.
                                                                    [required]
  --output-path ARTIFACT  Path where output artifact should be written.
                                                                    [required]
  --input-format TEXT     The format of the data to be imported. If not
                          provided, data must be in the format expected by the
                          semantic type provided via --type.
  --show-importable-types Show the semantic types that can be supplied to
                          --type to import data into an artifact.
  --show-importable-formats
                          Show formats that can be supplied to --input-format
                          to import data into an artifact.
  --help                  Show this message and exit.
```

クオリティチェック（オプション）

```
qiime demux summarize \
--i-data ${prefix}.qza \
--o-visualization ${prefix}.qzv
```

### 2) Denoising / OTU clustering

#### ASV

QIIME 2で利用できるノイズ除去の手法はDADA2とDeblur。DADA2にもDeblurキメラチェックとアバンダンスフィルタリングが含まれている。
- DADA2：エラー配列を真の配列に修正し、加算している。どのように修正するかをエラーモデルの学習によって推定している。
- Deblur：エラー配列を修正せずに除去する。DADA2 より配列の精度は高いが得られる read 数は少ない。先にクオリティスコアによるフィルタリングが必要。トリミング長は自動で設定される。

##### DADA2

DADA2のオプションを参照。

```
qiime dada2 denoise-paired \
--i-demultiplexed-seqs ${prefix}.qza \
--p-trunc-len-f 230 \ # DADA2のtruncLen。何サイクル目（何base目）までを使用するか。これより短いreadは自動的に捨てられます
--p-trunc-len-r 220 \
--p-trim-left-f 0 \ # 3'末をこのbase数トリムする（アダプター除去など？）
--p-trim-left-r 0 \
--p-trunc-q 2 \ # truncQ。スコアが指定した値以下の位置以降を切り捨て？
--p-max-ee-f 2.0 \ # maxEE。Trancationを行なった後に、指定した値より大きいexpected errorsを示した read を捨てます
--p-max-ee-r 2.0 \
--p-n-reads-learn 1000000 \ # エラーモデルの学習に要するread数
--p-n-threads 4 \
--o-table ${prefix}_table.qza \ # 各配列（ASV）の存在量データを出力する。塩基配列は含まれない。
--o-denoising-stats ${prefix}_denoiseStats.qza \ # DADA2 の各工程で入力された read 数を出力する。
--o-representative-sequences ${prefix}_repSeqs.qza \ # 各ASVの塩基配列を出力する。
--verbose
```

出力の.qzaを.qzvに変換すると[QIIME2 view](https://view.qiime2.org/)で見れるようになる

```
qiime metadata tabulate \
  --m-input-file    ${prefix}_denoiseStats.qza \
  --o-visualization ${prefix}_denoiseStats.qzv

qiime feature-table summarize \
  --i-table                ${prefix}_table.qza \
  --o-visualization        ${prefix}_table.qzv \
  --m-sample-metadata-file ${prefix}_metadata.tsv

qiime feature-table tabulate-seqs \
  --i-data          ${prefix}_repSeqs.qza \
  --o-visualization ${prefix}_repSeqs.qzv
```

```
qiime tools view ${prefix}_table.qza
```

##### Deblur

使っていないのでオプション設定はまだよくわからない

```
# ペアエンドのマージが含まれていないので、事前に`qiime vsearch`でマージしておく
qiime vsearch join-pairs   \
  --i-demultiplexed-seqs ${prefix}.qza \
  --o-joined-sequences ${prefix}_joined.qza
  
qiime deblur denoise-16S \
  --i-demultiplexed-seqs ${prefix}_joined.qza \
  --p-trim-length 250 \
  --p-mean-error 0.005 \
  --p-indel-prob 0.01 \
  --p-indel-max 3.0 \
  --p-jobs-to-start 36 \
  --o-table ${prefix}_table.qza \ # 各配列（ASV）の存在量データを出力する。塩基配列は含まれない
  --o-stats ${prefix}_denoiseStats.qza \ #  read 数を出力する
  --o-representative-sequences ${prefix}_repSeqs.qza \ # 各ASVの塩基配列を出力する
  --verbose
```


#### OTU

[Clustering sequences into OTUs using q2-vsearch](https://docs.qiime2.org/2022.8/tutorials/otu-clustering/)
[Qiime2でOTU法を行う](https://qiita.com/geed_rn/items/2bba7cb56a3f01dfc479)
[De novo VS Closed Ref](https://forum.qiime2.org/t/de-novo-vs-closed-ref/7024)

- de novo clustering
- Closed-reference Clustering
- Open-reference Clustering
の３種類がある。論文ではOpen-reference Clusteringが多い？

DADA2と違ってQCとペアマージが含まれないので別途行う必要がある。QIIME2内でも`qiime vsearch join-pairs` `qiime quality-filter q-score`でできるみたい。

- [Alternative methods of read-joining in QIIME 2](https://docs.qiime2.org/2020.11/tutorials/read-joining/)

キメラ除去について

- [Identifying and filtering chimeric feature sequences with q2-vsearch](https://docs.qiime2.org/2020.2/tutorials/chimera/)

Example
```
# ペアマージとマージ結果のsummarise
qiime vsearch join-pairs \
  --i-demultiplexed-seqs ${prefix}.qza \
  --o-joined-sequences   ${prefix}_joined.qza
qiime demux summarize \
  --i-data          ${prefix}_joined.qza \
  --o-visualization ${prefix}_joined.qzv
  
# クオリティフィルター
qiime quality-filter q-score \
  --p-min-quality 20 \
  --i-demux              ${prefix}_joined.qza \
  --o-filtered-sequences ${prefix}_joined_filtered.qza \
  --o-filter-stats       ${prefix}_joined_filter_stat.qza
qiime metadata tabulate \
  --m-input-file    ${prefix}_joined_filter_stat.qza \
  --o-visualization ${prefix}_joined_filter_stat.qzv
 
# FeatureTableとFeatureDataを作成
qiime vsearch dereplicate-sequences \
  --i-sequences              ${prefix}_joined_filtered.qza \
  --o-dereplicated-table     ${prefix}_table.qza \
  --o-dereplicated-sequences ${prefix}_repSeqs.qza

# OTU作成(open reference clustering)
qiime vsearch cluster-features-open-reference \
  --i-table               ${prefix}_table.qza \
  --i-sequences           ${prefix}_repSeqs.qza \
  --i-reference-sequences ！！ \
  --p-perc-identity 0.97 \
  --o-clustered-table         ${prefix}_table_or97.qza \
  --o-clustered-sequences     ${prefix}_repSeqs_or97.qza \
  --o-new-reference-sequences ${prefix}_newRepSeqs_or97.qza
  
# de novo キメラチェック、summarise
qiime vsearch uchime-denovo \
  --i-table     ${prefix}_table_pid97.qza \
  --i-sequences ${prefix}_repSeqs_pid975.qza \
  --output-dir  uchime_denovo
qiime metadata tabulate \
  --m-input-file    uchime_denovo/stats.qza \
  --o-visualization uchime_denovo/stats.qzv

# Exclude chimeras and “borderline chimeras”
qiime feature-table filter-features \
  --i-table ${prefix}_table_pid97.qza \
  --m-metadata-file  uchime_denovo/nonchimeras.qza \
  --o-filtered-table uchime_denovo/table-nonchimeric-wo-borderline.qza
qiime feature-table filter-seqs \
  --i-data ${prefix}_repSeqs_pid975.qza \
  --m-metadata-file uchime_denovo/nonchimeras.qza \
  --o-filtered-data uchime_denovo/rep-seqs-nonchimeric-wo-borderline.qza
qiime feature-table summarize \
  --i-table         uchime_denovo/table-nonchimeric-wo-borderline.qza \
  --o-visualization uchime_denovo/table-nonchimeric-wo-borderline.qzv
```


### 3) Taxonomy classification and taxonomic analyses

- [Example Workflow](https://github.com/otagoedna/getting_started_with_qiime2/blob/master/first_workflow.md)
- [Exploring_Taxonomy_Assignment](https://otagoedna.github.io/getting_started_with_qiime2/taxonomy_assignment/Exploring_Taxonomy_Assignment.html)
- [Overview of QIIME 2 Plugin Workflows - Taxonomy classification and taxonomic analyses](https://docs.qiime2.org/2022.8/tutorials/overview/#taxonomy-classification-and-taxonomic-analyses)

> classify-consensus-blast と classify-consensus-vsearch はアライメントに基づく分類法で、N 個のトップヒットの中からコンセンサスとなる分類を見つけます。これらの手法は、参照データベースである FeatureData[Taxonomy] と FeatureData[Sequence] ファイルを直接取得するため、事前学習は必要ありません。

> 機械学習ベースの分類法は classify-sklearn から利用でき、理論的には scikit-learn で利用可能な分類法のいずれかを適用することができます。これらの分類器は、例えば、各分類群を最もよく区別する特徴を学習するために学習する必要があり、分類プロセスに追加のステップが追加されます。分類器の学習は、参照データベースとマーカー遺伝子に固有であり、マーカー遺伝子と参照データベースの組み合わせごとに一度だけ行う必要があります。いくつかの学習済み分類器が提供されているので、通常は学習ステップを実行する必要はありません。

### Using Naive Bayes

```
qiime feature-classifier classify-sklearn \
  --i-classifier ~/db/qiime2/silva-138-99-nb-classifier.qza \
  --i-reads          ${prefix}_repSeqs.qza \
  --o-classification ${prefix}_taxonomyNB.qza

# visualization
qiime metadata tabulate \
  --m-input-file    ${prefix}_taxonomyNB.qza \
  --o-visualization ${prefix}_taxonomyNB.qzv
qiime taxa barplot \
  --i-table    ${prefix}_table.qza \
  --i-taxonomy ${prefix}_taxonomyNB.qza \
  --m-metadata-file ${prefix}_metadata.tsv \
  --o-visualization ${prefix}_taxonomyNB_barplot.qzv
```

#### Using BLAST

```
qiime feature-classifier classify-consensus-blast \
  --i-query.             ${prefix}_repSeqs.qza \
  --i-reference-reads    ~/db/qiime2/silva-138-99-seqs.qza \
  --i-reference-taxonomy ~/db/qiime2/silva-138-99-tax.qza \
  --p-perc-identity      0.97 \
  --o-classification     ${prefix}_taxonomyBlId97.qza \
  --verbose

# visualization (Naive Bayesと同じ)
qiime metadata tabulate \
  --m-input-file    ${prefix}_taxonomyBlId97.qza \
  --o-visualization ${prefix}_taxonomyBlId97.qzv
qiime taxa barplot \
  --i-table    ${prefix}_table.qza \
  --i-taxonomy ${prefix}_taxonomyBlId97.qza \
  --m-metadata-file ${prefix}_metadata.tsv \
  --o-visualization ${prefix}_taxonomyBlId97_barplot.qzv
```


#### Using vsearch

```
qiime feature-classifier classify-consensus-vsearch \
  --i-query              ${prefix}_repSeqs.qza \
  --i-reference-reads    ~/db/qiime2/silva-138-99-seqs.qza \
  --i-reference-taxonomy ~/db/qiime2/silva-138-99-tax.qza \
  --p-perc-identity 0.97 \
  --o-classification ${prefix}_taxonomyVsId97.qza \
  --verbose
```



## Example: Diversity Analyses

### 1) Sequence alignment and phylogeny building

#### 代表配列からtree作成する場合

[QIIME2docs - Phylogenetic inference with q2-phylogeny](https://docs.qiime2.org/2022.8/tutorials/phylogeny/)

`qiime phylogeny align-to-tree-mafft-fasttree`で以下４つのステップを実行できる。
- `qiime alignment mafft`：MAFFTを使ったアライメント
- `qiime alignment mask`：情報が少ない/曖昧なサイトをマスク
- `qiime phylogeny fasttree`：系統樹作成
- `qiime phylogeny midpoint-root`：系統樹を中点？でroot

```
# アライメント、マスク後アライメント、無根系統樹、有根系統樹を出力する
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences        ${prefix}_repSeqs.qza \
  --o-alignment        ${prefix}_repSeqs_aligned.qza \
  --o-masked-alignment ${prefix}_repSeqs_maskedAligned.qza \
  --o-tree             ${prefix}_unrootTree.qza \
  --o-rooted-tree      ${prefix}_rootedTree.qza

# または、出力ファイルのディレクトリ名だけを指定することもできる（上記４ファイルがディレクトリ内に出力される）
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences ${prefix}_repSeqs.qza \
  --output-dir  mafft_fasttree_output
```

#### 代表配列を既存treeへ配置する場合

短い16Sの場合、代表配列だけで系統樹を作るよりも、SEPP良いうのを使って既存系統樹の中に代表配列を配置する方が良いらしい。既存treeの配列と75%以上の総同性がない配列は捨てられる。
- [Ranking the biases: The choice of OTUs vs. ASVs in 16S rRNA amplicon data analysis has stronger effects on diversity measures than rarefaction and OTU identity threshold](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0264443#sec002)
- [Plugin website](https://github.com/qiime2/q2-fragment-insertion)
- [q2-fragment-insertion](https://library.qiime2.org/plugins/q2-fragment-insertion/16/)

[QIIME2 - Data resources](https://docs.qiime2.org/2022.8/data-resources/)からSEPP用レファレンスデータベースをダウンロード（系統樹？）

V1-V2領域human gutで見て見たところ、align-to-tree-mafft-fasttreeでの系統樹作成とfragment-insertionを比較したところ、unifrac距離の相関はweightedで0.892、unweightedで0.948だった。weighted Unifracは全体的に高くなっていた。Faith’s Phylogenetic Diversityも全体的に高くなっていて（当たり前？）、相関係数は0.988だった。

```
cd ~/db/qiime2/sepp
wget https://data.qiime2.org/2022.8/common/sepp-refs-silva-128.qza
wget https://data.qiime2.org/2022.8/common/sepp-refs-gg-13-8.qza
```

```
qiime fragment-insertion sepp \
  --p-threads 22 \
  --i-representative-sequences ${prefix}_repSeqs.qza \
  --i-reference-database ~/db/qiime2/sepp/sepp-refs-silva-128.qza \
  --o-tree ${prefix}_insertionTree.qza \
  --o-placements ${prefix}_insertionPlacements.qza
```

レファレンスtree内の配列と75%以上の相同性がない配列は捨てられるので、feature-tableもtreeに含まれた配列だけ残すようフィルターする。

```
qiime fragment-insertion filter-features \
  --i-table          ${prefix}_table.qza \
  --i-tree           ${prefix}_insertionTree.qza \
  --o-filtered-table ${prefix}_filteredTable.qza \
  --o-removed-table  ${prefix}_removedTable.qza \
  --verbose
  
# `--p-sampling-depth`を決める用
qiime feature-table summarize \
  --i-table                ${prefix}_filteredTable.qza \
  --o-visualization        ${prefix}_filteredTable.qzv \
  --m-sample-metadata-file ${prefix}_metadata.tsv
```

### 2) Core Diversity metricsの作成

`--p-sampling-depth`でrarefactionのsampling depthを設定する必要があるので、よき値を決めてください。各サンプルのリード数は、feature table（ここでは`${prefix}_table.qzv`）のtable.qzvのInteractive Sample Detail - Feature Countや、DADA2時のstats（ここでは`${prefix}_denoiseStats.qzv`）を見るとよい。指定リード数未満のサンプルは自動的に除かれるそうです。

作成されるファイルについて：[QIIME2 で生成するcore-metrics-resultsについて（１）](https://qiita.com/akari5/items/f25c16c5efaebe8427e4)


```
# mafft_fasttreeの系統樹を使う場合はこんな感じ
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny      ${prefix}_rootedTree.qza \
  --i-table          ${prefix}_table.qza \
  --p-sampling-depth ${integer} \
  --m-metadata-file  ${prefix}_metadata.tsv \
  --output-dir       ${prefix}_coreMetricsResult

# fragment_insertionのtreeを使う場合はこんな感じ
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny      ${prefix}_insertionTree.qza \
  --i-table          ${prefix}_filteredTable.qza \
  --p-sampling-depth 12000 \
  --m-metadata-file  ${prefix}_metadata.tsv \
  --output-dir       ${prefix}_coreMetricsResult_fragInsVer
```
 
${prefix}_coreMetricsResult内にたくさんファイルができる↓.qzvの図はemperorというなんだか見にくい感じのプロット。

Output artifacts(.qza):
- faith_pd_vector.qza: alpha-diversity。系統的な違い（系統樹の距離）から計算されるらしい。枝長の合計？
- unweighted_unifrac_distance_matrix.qza: beta-diversity。系統発生距離と在不在データから計算される。
- bray_curtis_pcoa_results.qza:
- shannon_vector.qza: alpha-diversity。種数(richness)を反映する。-Σ(pi*log(pi))（pはサンプル中のそのOTU(ASV)の割合）。QIIMEではlogの底が2なので注意（veganはe）。
- rarefied_table.qza:
- weighted_unifrac_distance_matrix.qza: beta-diversity。系統発生距離とカウントデータから計算される。
- jaccard_pcoa_results.qza:
- weighted_unifrac_pcoa_results.qza:
- observed_features_vector.qza: alpha-diversity。あるサンプルで観察されたOTU(ASV)の数。richness。
- jaccard_distance_matrix.qza: beta-diversity。生物の在不在データからコミュニティ間の多様性の違いを計算している。veganのjaccard(binaly=TRUE)と同じ。
- evenness_vector.qza: alpha-diversity。組成の均等度。Pielouのevenness。(shannon Index)/log(見つかったOTU数)。ちなみにlogをlogで割っているので底の違いの影響がなくなる。
- bray_curtis_distance_matrix.qza: beta-diversity。カウントデータデータを使った組成の非類似度。
- unweighted_unifrac_pcoa_results.qza:
Output visualizations(.qzv):
- unweighted_unifrac_emperor.qzv:
- jaccard_emperor.qzv:
- bray_curtis_emperor.qzv:
- weighted_unifrac_emperor.qzv:

ちなみに：Weighted UniFracは存在量を考慮するので、高頻度の分類群の影響を受けやすい。Unweighted UniFracは存在の有無の身を考慮して存在量を考慮しないので、低頻度の分類群の影響を受けやすい。[Chen et al. 2012](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3413390/)ではweited/unweited Unifrac距離はそれぞれ高頻度/低頻度分類群の影響を受けすぎるとしてgeneralized Unifrac距離を提唱している。[Unifracについて](https://github.com/biocore/unifrac)

### 3-a) Alpha diversity

`qiime diversity alpha-group-significance`: メタデータのカテゴリカルデータ列とAlpha diversityとの関連を検証する

```
# diversity
qiime diversity alpha-group-significance \
  --i-alpha-diversity ${prefix}_coreMetricsResult/faith_pd_vector.qza \
  --m-metadata-file   ${prefix}_metadata.tsv \
  --o-visualization   ${prefix}_coreMetricsResult/faith_pd_group_significance.qzv

# evenness
qiime diversity alpha-group-significance \
  --i-alpha-diversity ${prefix}_coreMetricsResult/evenness_vector.qza \
  --m-metadata-file   ${prefix}_metadata.tsv \
  --o-visualization   ${prefix}_coreMetricsResult/evenness_group_significance.qzv
```

`qiime diversity alpha-correlation`: メタデータの連続量列とAlpha diversityとの関連を検証する

```
qiime diversity alpha-correlation \
  --i-alpha-diversity ${prefix}_coreMetricsResult/faith_pd_vector.qza \
  --m-metadata-file   ${prefix}_metadata.tsv \
  --o-visualization   ${prefix}_coreMetricsResult/faith_pd_correlation_significance.qzv
```

### 3-b) Beta diversity

PERMANOVA。`-p-pairwise`でペアワイズ検定も行う。ここではメタデータのどの列を使うかも指定する必要がある。距離行列はどれでも（ここではunweighted_unifrac）。

```
qiime diversity beta-group-significance \
  --i-distance-matrix ${prefix}_coreMetricsResult/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file   ${prefix}_metadata.tsv \
  --m-metadata-column ${column_name} \
  --o-visualization   ${prefix}_coreMetricsResult/unweighted_unifrac_significance_${column_name}.qzv \
  --p-permutations 999
```


## Export data（qiime tools export）

.qzaからテキストファイルを出力する。コマンドは`qiime tools export`や`qiime tools extract`。ファイル名は指定できないみたい（`--output-path`で指定した値はディレクトリ名になり、その中に自動で名付けられたファイルが出力される）。

[QIIME2docs - Exporting data](https://docs.qiime2.org/2022.8/tutorials/exporting/)
[QIIME2 のデータをLEfSe 用に変更する](https://qiita.com/kohei-108/items/7bf2e59dd11b6fe4cbdb)


DADA2のstatsは.tsvになる

```
qiime tools export \
  --input-path    ${prefix}_denoiseStats.qza \
  --output-path   exported/denoiseStats # denoiseStatsフォルダ内にstats.tsvができる。
```

table(abundance)はbiom形式になるので、tsvにするにはその後textに変換する必要がある。

```
qiime tools export \
  --input-path ${prefix}_table.qza \
  --output-path exported/table # tableフォルダ内にfeature-table.biomができる
biom convert \
  -i exported/table/feature-table.biom \
  -o exported/table/feature-table.tsv \
  --to-tsv

>less exported/table/feature-table.tsv
# Constructed from biom file			
#OTU ID	u00005_201208	u00005_210107	u00005_210213
3016b070726d0a9e6cdb87c4a3914d07	0.0	0.0	0.0
44001a72292c6ac78f98d2ffaac23065	0.0	0.0	0.0
186d018791a835c7e871cf8c05caeedc	116.0	520.0	189.0
6af5f8ea781a98464482a500ff544d17	0.0	0.0	0.0
2bed2a2ff3e6c32e8c2a66e4832ac046	314.0	1230.0	981.0
021a9dd9ded3236d26bae49bdf2e7efd	0.0	0.0	0.0
5de4054ae7aa4f9985b43ac85e302867	0.0	0.0	0.0
3596a9491ac0346b7e2db1c16c03d8de	0.0	0.0	0.0
```


taxonomy情報。Feature ID、Taxon、Confidenceの3列のtsvとなる。

```
qiime tools export \
  --input-path    ${prefix}_taxonomyNB.qza \
  --output-path   exported/taxonomyNB # taxonomyNBフォルダ内にtaxonomy.tsvができる

>less exported/taxonomyNB/taxonomy.tsv
Feature ID      Taxon   Confidence
3016b070726d0a9e6cdb87c4a3914d07        d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Bacteroides; s__Bacteroides_dorei 0.995997222811799
44001a72292c6ac78f98d2ffaac23065        d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Bacteroides; s__Bacteroides_vulgatus      0.909845777736253
186d018791a835c7e871cf8c05caeedc        d__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Bifidobacteriales; f__Bifidobacteriaceae; g__Bifidobacterium; s__Bifidobacterium_pseudocatenulatum      0.9779750269072935
```

distance matrix。tsvになる（1列目の列名は空）。

```
qiime tools export \
  --input-path    ${prefix}_coreMetricsResult/weighted_unifrac_distance_matrix.qza \
  --output-path   exported/weighted_unifrac_distance # フォルダ内にdistance-matrix.tsvができる

>less exported/weighted_unifrac_distance/distance-matrix.tsv
	u00005_201208	u00005_210107	u00005_210213
u00005_201208	0	0.3333298735534015	0.2719313079038719
u00005_210107	0.3333298735534015	0	0.1948726923149328
u00005_210213	0.2719313079038719	0.1948726923149328	0
```

## Export data（qiime2R）

Rパッケージqiime2Rで.qzaをRに読み込めるので、そこからtsv出力もできる。read_qza()で読み込むとリストになる。以下主な要素。
- $uuidunique identifier for this object
- $type: the semantic type of the object
- $format: the format of the qiime artifact
- $data: the raw data ex OTU table as matrix or tree in phylo format

主な.qzaファイルを読み込んでデータ部分をtsvで出力

```
library(tidyverse)
library(qiime2R)

files <- c(
    "adults_table"                       = "~/mykinso/qiime2/adults_table.qza", # matrix
    "adults_taxonomyNB"                  = "~/mykinso/qiime2/adults_taxonomyNB.qza",
    "adults_coreMetricsResult/weighted_unifrac_distance_matrix"   = "~/mykinso/qiime2/adults_coreMetricsResult/weighted_unifrac_distance_matrix.qza", # distance matrix
    "adults_coreMetricsResult/unweighted_unifrac_distance_matrix" = "~/mykinso/qiime2/adults_coreMetricsResult/unweighted_unifrac_distance_matrix.qza", # distance matrix
    "adults_coreMetricsResult/shannon_vector"                     = "~/mykinso/qiime2/adults_coreMetricsResult/shannon_vector.qza"
)

data <- map(files, ~read_qza(.) %>% .$data) %>%
    set_names(names(files))

data[["adults_table"]] <- data[["adults_table"]] %>%
    as.data.frame() %>%
    rownames_to_column("ASV") %>%
    as_tibble()
data[["adults_coreMetricsResult/weighted_unifrac_distance_matrix"]] <- data[["adults_coreMetricsResult/weighted_unifrac_distance_matrix"]] %>%
    as.matrix() %>%
    as.data.frame()
data[["adults_coreMetricsResult/unweighted_unifrac_distance_matrix"]] <- data[["adults_coreMetricsResult/unweighted_unifrac_distance_matrix"]] %>%
    as.matrix() %>%
    as.data.frame()

for (i in 1:length(data)){
    name_i <- names(data)[i]
    data_i <- data[[i]]
    if (str_detect(name_i, "distance_matrix")){
        write.table(data_i, paste0("~/mykinso/qiime2/", name_i, ".tsv"), col.names = TRUE, row.names = TRUE, sep = "\t", quote=FALSE, na="")
    } else {
        write_tsv(data_i, paste0("~/mykinso/qiime2/", name_i, ".tsv"))
    }
}
```

## QIIME2結果をphyloseqに渡す

phyloseq.mdを参照。