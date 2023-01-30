# Rで系統樹を扱う

[biostatistics - 系統樹 ape ade4](https://stats.biopapyrus.jp/r/graph/phylogenetic-tree.html)

```
install.packages("ape")
```

パッケージを呼び出し

```
library(ape)
```

データ読み込み

```
file <- "~/mytreefile.tre"
tree <- read.tree(file)
```

中身はこんな感じ

```r
str(tree)

# List of 5
#  $ edge       : int [1:41, 1:2] 23 23 24 24 25 25 26 27 28 29 ...
#  $ edge.length: num [1:41] 0.00614 0.02231 0.03136 0.02861 0.02697 ...
#  $ Nnode      : int 20
#  $ tip.label  : chr [1:22] "W0004_v4|person-W0004|serial-v4" "W0004_v2|person-W0004|serial-v2" "W0004_v1|person-W0004|serial-v1" "W0005_v2|person-W0005|serial-v2" ...
#  $ root.edge  : num 0
#  - attr(*, "class")= chr "phylo"
#  - attr(*, "order")= chr "cladewise"
```

プロット

```
plot(
  tree,                    # read.treeで読み込んだデータ
  type = "phylogram",       # phylogram, cladogram, fan, unrooted, radial を指定することができる
  use.edge.length = TRUE,  # FALSE を指定すると、枝の長さは距離情報を含まなくなる！
  show.tip.label = TRUE,   # 葉のラベルを表示する
  show.node.label = TRUE,  # ノードのラベルを表示
  edge.color = "black",    # 枝の色
  edge.width = 1,          # 枝の太さ
  edge.lty = 1,            # 枝の種類（実線、点線など）
  root.edge = FLASE,       # 根を表示する
  tip.color = "black"      # 葉の色
)
```