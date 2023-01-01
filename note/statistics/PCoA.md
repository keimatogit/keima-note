# PCoA

## R

[バイオスタティスティクス基礎論 第 3 回 講義テキスト](https://lecture.ecc.u-tokyo.ac.jp/~aiwata/biostat_basic/2013/text4lec3.pdf)：cmdscale結果から寄与率の計算

distanceがeuclideanの時PCAと同じ。
でもprcompとはグラフの位置が対象（pc1, pc2の符号が逆？）になってたけど。
```r
suppressMessages(library(tidyverse))
suppressMessages(library(vegan))
s

mypcoa <- function(df, distance, suffix="_[0-9]+$", group=NA_character_, scale=FALSE, palette="Set2", shape=c(16,17,15,18,4,6,7,8), legend=TRUE){

    df <- as.data.frame(df)

    if(!is.numeric(pull(df,1))){
        rownames(df) <- df[,1]
    }

    df_t <- df %>%
        select(where(is.numeric)) %>%
        t() %>%
        as.data.frame()

    if(all(is.na(group))){
        group <- str_remove(rownames(df_t), suffix)
    } else {
        stopifnot("'group' must be the same length as samples (columns containing numeric values)." = (length(group) == nrow(df_t)))
    }

    dist_mat <- vegdist(df_t, method=distance)
    pcoa_result <- cmdscale(dist_mat, k= 5, eig=TRUE) # PCo5まで

    # スコア
    PCo_score <- pcoa_result$points %>%
        as_tibble() %>%
        rename(1:5 %>% set_names(paste0("PCo", 1:5))) %>%
        mutate(sample = rownames(df_t),
               group  = group) %>%
        select(sample, group, everything())

    # 寄与率
    contribution_ratio <- pcoa_result$eig[1:5] / sum(pcoa_result$eig)
    names(contribution_ratio) <- paste0("PCo", 1:5)

    # グラフ

    xlabel <- paste0("PCo1 (", (contribution_ratio["PCo1"]*100) %>% round(1), "%)")
    ylabel <- paste0("PCo2 (", (contribution_ratio["PCo2"]*100) %>% round(1), "%)")
 
    g <- ggplot(PCo_score, aes(x=PCo1, y=PCo2, color=group, shape=group)) +
        scale_color_brewer(palette = palette) +
        scale_shape_manual(values = shape) +
        labs(x = xlabel, y = ylabel) +
        theme(text = element_text(size = 18),
              legend.title = element_blank(),
              panel.background = element_rect(fill = "white", colour = "black")) +
        coord_fixed()

    if(legend == FALSE){
        g <- g + guides(color = "none", shape = "none")
    }

    g1 <- g +
        geom_point(size=2.4)

    g2 <- g +
        geom_text(aes(label = sample))

    return(
        list(
            result    = pcoa_result,
            PCo_score = PCo_score,
            contribution_ratio = contribution_ratio,
            plots = list(group = g1, sample = g2)
        )
    )

}

```