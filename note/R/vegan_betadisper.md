# vegan - betadisper

PERMDISP：グループ間の分散の違いを検定。PERMANOVA（平均（重心）のほかに分散の違いにも影響される）と併用することが多い。

[betadisper examples for EqCov paper](https://mattsigal.github.io/eqcov_supp/betadisp-ex.html)

引数はdistance matrix（データフレームをvegdistで取得）。
example: `betadisp_result <- betadisper(distance_matrix, group)`
permutationテストは、`permutest(betadisp_result)`

```
library(vegan)

group <- colnames(data_g)[-1] %>%
  str_remove("[0-9]{2}") %>%
  factor(., levels = c("V","C"))

dist_methods <- c("euclidean", "bray")

distmat_list <- map(dist_methods, function(dist_i){
  data_g_cpm %>%
   filter(taxname != "unannotated_Bacteria") %>%
   column_to_rownames("taxname") %>%
   t() %>%
   vegdist(method = dist_i)
}) %>% set_names(dist_methods)

disp_res_lst <- distmat_list %>% map(betadisper, group=group, type="median")
disp_res_permutest_lst <- map(disp_res_lst, ~permutest(., permutations = 10000))
```

プロット

```
plot(disp_res_lst[[1]], ellipse=TRUE, hull=TRUE)
```


寄与率入りプロット
参考：[betadisper examples for EqCov paper](https://mattsigal.github.io/eqcov_supp/betadisp-ex.html)

```
myplot <- function(disp_res, title){
    iris.bd <- disp_res

    # https://mattsigal.github.io/eqcov_supp/betadisp-ex.html　からコピペ
    labs <- paste("Dimension", 1:4, "(", 
                round(100*iris.bd$eig / sum(iris.bd$eig), 2), "%)")
    plot(iris.bd, cex=1, pch=15:17,
      main=title, cex.lab=1.25,
      xlab=labs[1], ylab=labs[2],
      hull=FALSE, ellipse=TRUE, conf=0.68, lwd=2)
}

options(repr.plot.width=10, repr.plot.height=14)
split.screen(c(2,1)) # 2行1列
screen(1)
myplot(disp_res_lst[[1]], "Euclidean")
screen(2)
myplot(disp_res_lst[[2]], "Bray-Curtis")
```