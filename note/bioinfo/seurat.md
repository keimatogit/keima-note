---
---

# Seurat

> Seuratは、シングルセルのRNA-seqデータのQC、解析、および探索のために設計されたRパッケージです。

参考サイト
- [Seurat](https://satijalab.org/seurat/)
- [Seurat - get started](https://satijalab.org/seurat/articles/get_started.html)
- [Seurat - Guided Clustering Tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)
- [Seurat - reference](https://satijalab.org/seurat/reference/index.html)

略語
- UMI: Unique moleqular identifier
- HTOs: hashtag oligos
- ADT: antibody-derived tags

## インストール

condaで

```
conda install -c bioconda r-seurat
```

Rで

```
install.packages("Seurat")
```

## Example

### パッケージ読み込み

```
library(tidyverse)
library(Seurat)
library(RColorBrewer)
```

### Reading data

[Load in data from 10X](https://satijalab.org/seurat/reference/read10x)

cellrangerの出力フォルダを指定してデータを読み込みます。Gene Expressionには通常のscRNA発現量データ、Antibody Captureには表面タンパクの発現量データが格納されています。

```
data_dir <- "/imetgpfs/projects/chip/v562/F4199/TIL_WT1-CTL_5DE/outs/filtered_feature_bc_matrix/"
data <- Read10X(data.dir = data_dir)

# 中身はこんな感じ
str(data)
# List of 2
#  $ Gene Expression :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#   .. ..@ i       : int [1:6732305] 10 16 17 20 38 40 52 69 96 131 ...
#   .. ..@ p       : int [1:2758] 0 2655 5933 6114 8029 9245 13434 13529 13655 14673 ...
#   .. ..@ Dim     : int [1:2] 32285 2757
#   .. ..@ Dimnames:List of 2
#   .. .. ..$ : chr [1:32285] "Xkr4" "Gm1992" "Gm19938" "Gm37381" ...
#   .. .. ..$ : chr [1:2757] "AAACCTGCAATCACAC-1" "AAACCTGCACCAGATT-1" "AAACCTGCAGGACCCT-1" "AAACCTGCAGTCGTGC-1" ...
#   .. ..@ x       : num [1:6732305] 4 1 1 1 1 1 1 1 1 2 ...
#   .. ..@ factors : list()
#  $ Antibody Capture:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#   .. ..@ i       : int [1:8269] 0 1 2 0 1 2 0 1 2 0 ...
#   .. ..@ p       : int [1:2758] 0 3 6 9 12 15 18 21 24 27 ...
#   .. ..@ Dim     : int [1:2] 3 2757
#   .. ..@ Dimnames:List of 2
#   .. .. ..$ : chr [1:3] "Hashtag1" "Hashtag2" "Hashtag3"
#   .. .. ..$ : chr [1:2757] "AAACCTGCAATCACAC-1" "AAACCTGCACCAGATT-1" "AAACCTGCAGGACCCT-1" "AAACCTGCAGTCGTGC-1" ...
#   .. ..@ x       : num [1:8269] 25 7744 50 19 84 ...
#   .. ..@ factors : list()
```

seurat objectを作成します。seurat objectはS4クラスで、例えばseurat_object@assaysに発現量データ、seurat_object@meta.dataに色んなメタデータが入ります。またseurat_object@active.assayはデフォに設定したアッセイで、色んな関数を実行する際にアッセイ名を省略した場合はこのアッセイに適応されます（通常はRNA）。

メモ１: seurat_object@assays$RNAは、単に`seurat_object[["RNA"]]`や`seurat_object$RNA`で呼び出し可能。

メモ２: seurat_object@meta.dataはデータフレーム。各列の呼び出しは、例えばnCount_RNAの場合には、単に`seurat_object$nCount_RNA`（named vector）や`seurat_object[["nCount_RNA"]]`（データフレームの列抽出）で可能。

```
# Gene Expression
seurat_object = CreateSeuratObject(counts = data$`Gene Expression`)

# Antibody Capture
# Antibody Captureの方の名前（アッセイ名）は、HTO、ADTなどその時のデータに沿う名前をつけておきます
seurat_object[['HTO']] = CreateAssayObject(counts = data$`Antibody Capture`)

# 中身を確認しましょう
str(seurat_object)
```


### Filtering cells

質のよくない細胞データを除きます。チュートリアルに沿って、見つかった遺伝子数が少ないセルと、ミトコンドリア遺伝子割合が多いセルを除くことにします。

まずはミトコンドリア遺伝子のパーセンテージを取得します。[PercentageFeatureSet()](https://satijalab.org/seurat/reference/percentagefeatureset)では、単純に特定のパターンが含まれる遺伝子名（feature名）を探して、それらのカウントの和を全ての遺伝子のカウント総和で割って100を掛けます。ミトコンドリア遺伝子を探す場合、パターンがmt-だったりMT-だったりするようなので注意（遺伝子名一覧は`rownames(seurat_object[['RNA']])`で見れるので、`sum(str_detect(rownames(seurat_object[['RNA']]), "^mt-"))`などでパターンが含まれる遺伝子名があるか確認しよう）。

```
# seurat_object@meta.data$percent.mtに結果が格納されます
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^mt-")
```

バイオリンプロットで、細胞ごとのfeature数（見つかった遺伝子数）、総カウント数、ミトコンドリア遺伝子割合をチェック。featuresで指定しているのは、seurat_object@meta.dataの列名です。

```
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

カウント数 x ミトコン遺伝子%と、カウント数 x Feature数の散布図

```
plot1 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
```

バイオリンプロットや散布図を参考に閾値を決めて、基準を満たさない細胞を除きます（基準を満たす細胞だけを抽出します）。

```
# feature数が1001から5999、かつミトコン遺伝子が5%未満の細胞だけを抽出
cutoff <- c("small" = 1000, "large"=6000, "mt"=5)
seurat_object <- subset(
    seurat_object,
    subset = nFeature_RNA > cutoff["small"] & nFeature_RNA < cutoff["large"] & percent.mt < cutoff["mt"]
)
```

フィルタリング後のバイオリンプロットをチェック。

```
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```


### Demultiplexing with hashtag

[Demultiplexing with hashtag oligos (HTOs)](https://satijalab.org/seurat/articles/hashing_vignette.html)

複数のサンプルをhashtagと呼ばれる抗体で染色して、混合して1サンプルとしてシーケンスした場合、hashtagに基づいて分離します（Cell Hashing）。どのハッシュタグの発現も不十分だったり（ネガティブ）、十分な発現が見られたハッシュタグが複数あった（ダブレット）細胞は除き、ひとつのハッシュだけ十分発現していた（シングレット）細胞だけ抽出します。また、どのハッシュタグに分類されるかの情報を追加します。

データ読み込み時、表面タンパクの発現データは"HTO"に入れました。ハッシュタグ名を確認しましょう。

```
rownames(seurat_object[["HTO"]])
```

発現データのノーマライズ（CLR: centered log-ratio transformation）

```
seurat_object <- NormalizeData(seurat_object, assay = "HTO", normalization.method = "CLR")
```

demultiplex。[HTODemux()](https://satijalab.org/seurat/reference/htodemux)では、まずそれぞれのハッシュタグについて陽性/陰性の区別をする閾値を計算して、いずれのハッシュタグでも閾値を超えなかった場合はNegative、複数ハッシュタグで閾値を超えていた場合にdoubletに分類されるそうです。（positive.quantileは閾値を決めるときの基準で、デフォは0.99です。これを下げると陽性が増えるのでNegativeは減りますがDoubletは増えます。）

```
seurat_object <- HTODemux(seurat_object, assay = "HTO", positive.quantile = 0.99)
```

結果を確認しましょう。seurat_object@meta.dataに、HTO_maxID、HTO_secondID、HTO_margin、HTO_classification、HTO_classification.global、hash.IDという列が追加されています。よく見るのは以下３つ。
- HTO_classification.global: Singlet or Doublet or Negativeの判定
- hash.ID: ハッシュタグ名 or Doublet or Negative
- HTO_classification: ハッシュタグ名 or Negative（DoubletはHashtag1_Hashtag2のように検出された複数タグ名が連結されて表示）

``` 
# 分類結果（細胞数）を確認
table(seurat_object$HTO_classification.global)
table(seurat_object$hash.ID)
table(seurat_object$HTO_classification)
```

各ハッシュタグの発現量の、HTO_maxID（たぶん一番発現が多かったタグ名）ごとのリッジプロット。ややこしいけど、Identity(Y軸)が細胞のグループ（Idents(seurat_object)で、ここではHTO_maxID）で、タイトルがハッシュタグ名(featuresで指定)。

```
num_of_HTO <- nrow(seurat_object[["HTO"]]@data) # ハッシュタグが何種類あるか
Idents(seurat_object) <- "HTO_maxID"

RidgePlot(seurat_object, assay = "HTO", features = rownames(seurat_object[["HTO"]])[1:num_of_HTO], ncol = num_of_HTO)
```

（補足）細胞グループごとに密度プロットを見たかったので自前プロット。HTO_classificationごとに各ハッシュタグ発現量の分布をプロット。

```
tmp <- seurat_object[["HTO"]]@data %>%
    t() %>%
    as_tibble(rownames="cell") %>%
    mutate(HTO_classification = seurat_object$HTO_classification) %>%
    pivot_longer(-c(cell, HTO_classification), names_to="hashtag", values_to="value")

g <- ggplot(tmp) +
    geom_density(aes(value, color=hashtag), size= 1) +
    scale_color_manual(values=cols$hashtag %>% unname()) +
    facet_wrap(. ~ HTO_classification, ncol=4, scales="free_y") +
    theme(text = element_text(size = 18),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = NULL),
        axis.line = element_line(lineend = "square"))

plot(g)
```

Doublet, Negative, Singletに分類された細胞のRNAカウント数の分布をバイオリンプロットでチェック。

```
Idents(seurat_object) <- "HTO_classification.global"
options(repr.plot.width=7, repr.plot.height=7)
VlnPlot(seurat_object, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
```

Doublet, Negativeを除いてsingletだけ抽出する

```
Idents(seurat_object) <- "HTO_classification.global"
seurat_object <- subset(seurat_object, idents = "Singlet")
# または
seurat_object <- subset(seurat_object, subset = HTO_classification.global == "Singlet")
```

DoubletとNegativeがなくなったことを確認しましょう

```
table(seurat_object$hash.ID)
```


### Normalizing data

RNAデータのノーマライズ

アッセイ名を指定しない場合に適応されるデフォルトアッセイ名を確認しておきます。ここでは'RNA'と'HTO'を作成していましたが、デフォルトアッセイは'RNA'のはず。
```
DefaultAssay(seurat_object)
```

ノーマライズ

```
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000) 
```

### Feature selection

細胞間での発現変動が大きい遺伝子を抽出します（デフォルトでは上位2000遺伝子を抽出）。これらの遺伝子セットに注目してPCAなど下流の解析が行われます。選ばれた遺伝子名はseurat_object[["RNA"]]@var.featuresに格納されます。VariableFeatures(seurat_object)で呼び出せます。

```
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
```

トップ10をチェック

```
top10 <- head(VariableFeatures(seurat_object), 10)
top10
```

平均発現量 x 分散(?)の散布図（トップ10をラベル）

```
plot1 <- VariableFeaturePlot(seurat_object)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
options(repr.plot.width=14, repr.plot.height=7)
plot1 + plot2
```


### Scaling

PCAなどのために、発現量の平均値がゼロ、細胞間の分散が1になるように発現量を変換します（高発現の遺伝子の影響が優位にならないようにするため）。変換後のデータはseurat_object[["RNA"]]@scale.dataに格納されます。

```
all.genes <- rownames(seurat_object)
seurat_object <- ScaleData(seurat_object, features = all.genes)
```

### PCA

PCAの実行。デフォルトではFeature selectionで抽出した遺伝子を使用します。npcsのデフォは50（PC50まで出す）。

```
seurat_object <- RunPCA(seurat_object, npcs=50)
```

PC1, PC2の各遺伝子の因子負荷量（デフォでは大きい順に30遺伝子分）をプロット。

```
VizDimLoadings(seurat_object, dims = 1:2, reduction = "pca")
```

PCAのプロット

```
DimPlot(seurat_object, reduction = "pca")
```


### Determine the ‘dimensionality’

細胞のクラスタリングはPCAのスコアを用いて行われるので、第何PCまで使用するかを決めます。

> ユーザーの皆様には、異なる数のPC（10、15、あるいは50！）でダウンストリーム解析を繰り返すことをお勧めします。その結果、多くの場合、結果は劇的に変わることはありません。
> このパラメータを選択する際は、高いほうを選ぶことをお勧めします。例えば、5個のPCでダウンストリーム解析を行った場合、結果に大きな悪影響を及ぼします。

「PCA の有意性についての JackStraw 解析の結果をプロット」（データ量が多い場合は時間がかかるので下のElbowPlot()が良いそうです）。よく分からないけど、点線が一様分布（null分布）で、そこに近い成分はもうあまり「有意ではない」そうです。それまでのPCと比べて急に下にきている線があればそこで切っても良い？

```
seurat_object <- JackStraw(seurat_object, num.replicate = 100)
seurat_object <- ScoreJackStraw(seurat_object, dims = 1:20)
JackStrawPlot(seurat_object, dims = 1:15)
```

「エルボープロット」：主成分の標準偏差をプロット。各々によって説明される分散の割合？がだんだん下がるので、平らになった後の部分（エルボー）はいらないみたい。

```
ElbowPlot(seurat_object, ndims=50)
```

PCいくつまで使用するかをnumPCに格納しておきます

```
numPC <- 20
```


### Determine the number of clusters

クラスタリングのパラメータ（resolution）を変えた時に、クラスター数がどう変わるかや、細胞がどのようにクラスター間を移動するかを可視化します。クラスタが互いにどのように関連しているか、どれが明確に区別され、どれが不安定であるかを見て、クラスタリングのクラスター数やパラメータを決定します。ここではseuratのほかにclustreeというパッケージを使います。

> resolution: 値が大きくなるほどクラスタの数が多くなります。このパラメータを0.4-1.2の間に設定すると、3K細胞程度の単一細胞データセットでは、一般的に良い結果が得られることがわかりました。より大きなデータセットでは、最適な解像度はより高くなることが多い。クラスターはIdents()関数を使って見つけることができます。

[Plotting clustering trees](https://cran.r-project.org/web/packages/clustree/vignettes/clustree.html)

```
library(clustree)

# Select a range of resolutions
resolution.range <- seq(from = 0.4, to = 1.2, by = 0.1)

# Find clusters using a range of resolutions
clustree_object <- seurat_object
clustree_object <- FindNeighbors(clustree_object, dims = 1:numPC)
clustree_object <- Seurat::FindClusters(clustree_object, resolution = resolution.range)

# 結果の表示
clustree(clustree_object, prefix = "RNA_snn_res.")
```


### Clustering the cells

決めたresolutionでクラスタリングします（resolutionデフォルトは0.8、チュートリアルでは0.5でした）。クラスタリング結果はmeta.data$RNA_snn_res.[resolution]に格納され、さらにmeta.data$seurat_clustersと、Idents(seurat_object)にも格納されます。meta.data$seurat_clustersとIdents(seurat_object)はクラスタリングをするたびに書き換えられるので、resolutionを変えて何度も実行した場合などは、最新の結果に置き換わります。

```
resolution <- 0.6
seurat_object <- FindNeighbors(seurat_object, dims = 1:numPC)
seurat_object <- FindClusters(seurat_object, resolution = resolution)

# 最初の5細胞のクラスターIDをチェック
head(Idents(seurat_object), 5)

# 各クラスターの細胞数をチェック
table(Idents(seurat_object))
```


### UMAP

`DimPlot()`ではlabel = TRUEでクラスター名をグラフ上に載せられる。プロットの出力はggplotなので、色のカスタムなど追加できる。

```
seurat_object <- RunUMAP(seurat_object, dims = 1:numPC)

# クラスターごとに色分けしたUMAP
umap1 <- DimPlot(seurat_object, reduction = "umap")

# ハッシュタグごとに色分けしたUMAP
umap2 <- DimPlot(seurat_object, reduction = "umap", group.by = "HTO_classification")

# 色のカスタム例
cols <- list(
    "cluster" = c(brewer.pal(9, "Oranges")[c(3,5,8)], brewer.pal(9, "Blues")[c(3,5,8)], brewer.pal(9, "RdPu")[5]) %>%
        set_names(paste0("Cluster", c(1,2,5,0,3,6,4))),
    "hashtag" = brewer.pal(3, "Set2") %>% set_names(paste0("MyHashtag", 1:3))
)
umap1 <- umap1 +
    scale_color_manual(
        values = cols$cluster %>% set_names(names(cols$cluster) %>% str_remove("Cluster") %>% as.numeric),
        labels = names(cols$cluster)
    )
umap2 <- umap2 +
    scale_color_manual(
        values = cols$hashtag %>% set_names(paste0("Hashtag", 1:3))
        labels = names(cols$hashtag)
    )
```


### （補足）RDataの保存

Seuratの読み込みやいろんな関数の実行は時間がかかるので、このあたりででもseurat_objectをRDataで保存しておくと続きが読み出せて便利。

```
saveRDS(seurat_object, file = "seurat_object.rds")
```

呼び出す時

```
rds <- "seurat_object.rds"
seurat_object <- readRDS(rds)
```

### Finding DEFs

feature（遺伝子）発現量をクラスター間で比較して、発現量が異なる遺伝子を探します。デフォルトでは、比較方法はノンパラのWilcoxon Rank Sum testです。

[FindMarkers()](https://satijalab.org/seurat/reference/findmarkers)の主なオプション
- min.pct:どちらかのグループでmin.pct(0-1)以上の細胞で発現が見られる遺伝子のみ比較に使用します。デフォは0.1。
- logfc.threshold: 平均log2FoldChangeがlogfc.threshold以上の遺伝子のみ比較に使用します。デフォは0.25。

[Finding differentially expressed features (cluster biomarkers)](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#finding-differentially-expressed-features-cluster-biomarkers-)

```
# クラスター0とその他の細胞を比較
cluster0_markers <- FindMarkers(seurat_object, ident.1 = 0)

# クラスター1,2とその他の細胞を比較
cluster12_markers <- FindMarkers(seurat_object, ident.1 = c(0,1,2)

# クラスター0とクラスター1,2を比較
cluster0_12_markers <- FindMarkers(seurat_object, ident.1 = 0, ident.2 = c(1,2))
```

結果はデータフレームで、行名が遺伝子名、列は比較方法に応じた統計量で、wilcoxの場合は以下の５つ。出力時は特にソートされていないようなので、dplyr::arrangeでFoldCahnge降順などに並べ替えると良い。
- p_val: p値
- avg_log2FC: 平均log2FoldChange
- pct.1: 第1グループでその遺伝子が検出された細胞の割合
- pct.2: 第2グループでその遺伝子が検出された細胞の割合
- p_val_adj: 多重比較補正したp値（Bonferroni。「事前に遺伝子をフィルタリングしてテスト数を減らすため他の補正方法は推奨されません」）

[FindAllMarkers()](https://satijalab.org/seurat/reference/findallmarkers)は、各クラスターについて、そのクラスターとその他の細胞の比較を全部行います。主なオプションはFindMarkers()と同じ。

```
cluster_markers_all <- FindAllMarkers(seurat_object)
```

結果はデータフレームで、FindMarkers()の列に加えてclusterという列ができます。また、遺伝子名が行名だけでなくgeneという列に入ります。


### 遺伝子発現プロット

（ちなみにオプション`cols`で色をカスタムできます）

特定の遺伝子（ここではCd69, Lag3, Itgae）について、クラスターごとの発現量のバイオリンプロット

```
VlnPlot(seurat_object, features = c("Cd69", "Lag3", "Itgae"))

# 色カスタム例
cols <- list(
    "cluster" = c(brewer.pal(9, "Oranges")[c(3,5,8)], brewer.pal(9, "Blues")[c(3,5,8)], brewer.pal(9, "RdPu")[5]) %>%
        set_names(paste0("Cluster", c(1,2,5,0,3,6,4)))
)
VlnPlot(seurat_object, features = c("Cd69", "Lag3", "Itgae"), cols = cols$cluster %>% set_names(names(cols$cluster) %>% str_remove("Cluster") %>% as.numeric))
```

特定の遺伝子の発現量をUMAP上にカラーリングしたもの

```
FeaturePlot(seurat_object, features = c("Cd69", "Lag3", "Itgae"))
```


### （補足）クラスターxハッシュタグの棒グラフ

各クラスターの細胞には各ハッシュタグがどれくらいの割合で含まれているか？と、各ハッシュタグの細胞には各クラスターがどのくらいの割合で含まれているか？の積み上げ棒グラフ。

```
# データ整形

cluster_order <- c(1,2,5,0,3,6,4)

data_g <- tibble(
    hashtag = seurat_object@meta.data$hash.ID,
    cluster = seurat_object@meta.data$seurat_clusters
)
cluster_num <- length(unique(data_g$cluster))
data_g <- data_g %>%
    mutate(cluster = paste0("Cluster", cluster)) %>%
    mutate(hashtag = factor(hashtag, levels = paste0("Hashtag", 1:3), labels=paste0("C030", 1:3))) %>%
    mutate(cluster = factor(cluster, levels=paste0("Cluster", cluster_order)))
head(data_g)
```

```
# 棒グラフ用関数
my_barplot <- function(data, x, fill, position="fill"){
    if(position == "fill"){
        ylabel = "ratio"
    } else {
        ylabel = "count"
    }
    data <- data %>%
        mutate(!!fill := factor(get(fill), levels = rev(levels(get(fill)))))
    g <- ggplot(data) +
        geom_bar(aes(x = get(x), fill = get(fill)), position = position) +
        coord_flip() +
        scale_fill_manual(values = cols[[fill]]) +
        scale_x_discrete(limits = rev(levels(pull(data, !!x)))) +
        scale_y_continuous(expand = c(0, 0)) +
        labs(x = NULL, y = ylabel) +
        # guides(fill = guide_legend(reverse = TRUE)) +
        theme(text = element_text(size = 18),
            legend.title = element_blank(),
            panel.background = element_rect(fill = "white", colour = NULL),
             axis.line = element_line(lineend = "square"))
    return(g)
}

# カラーセットの準備
cols <- list(
    "cluster" = c(brewer.pal(9, "Oranges")[c(3,5,8)], brewer.pal(9, "Blues")[c(3,5,8)], brewer.pal(9, "RdPu")[5]) %>%
        set_names(paste0("Cluster", c(1,2,5,0,3,6,4))),
    "hashtag" = brewer.pal(3, "Set2") %>% set_names(paste0("C030", 1:3))
)
```

```
# 各クラスター細胞に各ハッシュタグがどのくらい含まれているか（g1: ratio, g2: count）
g1 <- my_barplot(data_g, x="cluster", fill="hashtag")
g2 <- my_barplot(data_g, x="cluster", fill="hashtag", position="stack")
g1 + g2
```

```
# 各ハッシュタグ細胞に各クラスターがどのくらい含まれているか（g1: ratio, g2: count）
g1 <- my_barplot(data_g, x="hashtag", fill="cluster")
g2 <- my_barplot(data_g, x="hashtag", fill="cluster", position="stack")
g1 + g2
```