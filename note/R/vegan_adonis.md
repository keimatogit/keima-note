# vegan - adonis

PERMANOVA：構成の平均（重心）の違いの検定（分散の違いにも影響される！）

引数はデータフレームまたはdistnace matrix


Example

```
dist_methods <- c("euclidean", "bray")

distmat_list <- map(dist_methods, function(dist_i){
  data_g_cpm %>%
   filter(taxname != "unannotated_Bacteria") %>%
   column_to_rownames("taxname") %>%
   t() %>%
   vegdist(method = dist_i)
}) %>% set_names(dist_methods)

group <- colnames(data_g)[-1] %>%
  str_remove("[0-9]{2}") %>%
  factor(., levels = c("V","C"))

anova_res_lst <- map(distmat_list, function(distmat_i){    
    adonis(formula = distmat_i ~ group, permutations = 10000)
}) 
```