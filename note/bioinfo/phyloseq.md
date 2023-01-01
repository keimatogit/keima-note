# phyloseq

[phyloseq](https://joey711.github.io/phyloseq/)

QIIME2とphyloseqのUnifrac距離は随分違うらしい
- [very different weighted unifrac values for qiime2 versus phyloseq](https://github.com/joey711/phyloseq/issues/956)


```
conda install -c bioconda bioconductor-phyloseq
```

#

phyloseqはS4 class。スロットは以下
- @otu_table; otu_table();   OTS/ASVカウントテーブル。matrix?
- @sam_data;  sample_data(); サンプルメタデータ。data.frame?
- @tax_table; tax_table();   分類アノテーション。matrix?
- @phy_tree;  phy_tree();    系統樹。
- @refseq;    refseq();      代表配列データ。DADA2のチュートリアルではBiostrings::DNAStringSetでS4クラスDNAStringSet


## DADA2の結果を使う

DADA2結果からphyloseqオブジェクトを作成（）。

```
# 
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), # サンプル x OTU(ASV)のmatrix
               sample_data(samdf), # 行名がサンプル名、各列がサンプル情報（任意）のデータフレーム
               tax_table(taxa)) # OTU(ASV) x taxonomy("Kingdom" "Phylum", ..., "Species")のmatrix
```

## QIIMEの結果を使う

[phyloseq (+α) で16S ampliconデータから色々な図を作ろう](https://qiita.com/xvtyzn/items/d7b0e6e2951f59e73f31)

phyloseqオブジェクトを作成（）。

```
library(tidyverse)
library(qiime2R)

workdir <- "/user/ngsdata/kmatsumoto/mykinso/qiime2/"

tree     <- read_qza(paste0(workdir, "adults_rootedTree.qza")) #rootつき系統樹の読み込み
asv_tb   <- read_qza(paste0(workdir, "adults_table.qza")) #カウントテーブルの読み込み
taxonomy <- read_qza(paste0(workdir, "adults_taxonomyNB.qza")) #分類アノテーションの読み込み
repseq   <- read_qza(paste0(workdir, "adults_repSeqs.qza")) #代表配列の読み込み
metadata <- read_tsv(paste0(workdir, "adults_metadata.tsv"), show_col_types = FALSE) %>% as.data.frame() %>% column_to_rownames("sample-id") #メタデータの読み込み

tax_table <- parse_taxonomy(taxonomy$data) # これもqiime2Rの関数
tax_table <- as.matrix(tax_table)

metadata <- sample_data(metadata) # これはphyloseqの関数

myphyloseq <- phyloseq(
    otus = otu_table(asv_tb$data, taxa_are_rows = T),
    tree = tree$data,
    tax_table = tax_table(tax_table),
    metadata = metadata,
    refseq = repseq$data
)
```

tree描画（代表100配列）

```
subGP1 <- prune_taxa(taxa_names(myphyloseq)[1:100], myphyloseq)
pdf("~/mykinso/qiime2/phyloseq_plot.pdf")
plot_tree(subGP1, color="user_id", justify= "justify")
dev.off()
```

unifrac計算 

```
# rarefaction 
myphyloseq_rarefied <- myphyloseq %>% 
    rarefy_even_depth(rngseed=1, sample.size=0.9*min(sample_sums(myphyloseq)), replace=F)

# weighted（normalizedはTRUE/FALSE）
phyloseq_weighted <- phyloseq::UniFrac(myphyloseq_rarefied, weighted=TRUE, normalized = FALSE)

# unweighted（normalizedはTRUE/FALSE）
phyloseq_unweighted <- phyloseq::UniFrac(myphyloseq_rarefied, weighted=FALSE, normalized = TRUE)
```