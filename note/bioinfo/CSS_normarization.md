# Cumulative Sum Scaling (CSS)

- [CSS - Metagenomics wiki](https://www.metagenomics.wiki/tools/16s/norm/css)
- [manual(pdf)](https://bioconductor.org/packages/release/bioc/manuals/metagenomeSeq/man/metagenomeSeq.pdf)


```
conda install -c bioconda megahit
```

使い方
```
# convert table into package format
metaSeqObject <- newMRexperiment(df)

# CSS normalization
metaSeqObject_CSS <- cumNorm(
  metaSeqObject,
  p = cumNormStatFast(metaSeqObject)
)

# convert CSS normalized data into data.frame-formatted OTU table (log2 transformed data)
df_CSS <- data.frame(MRcounts(metaSeqObject_CSS, norm=TRUE, log=TRUE))

```

- `metaSeqObject_CSS`では`Default value being used.`と表示される。`p = cumNormStatFast(metaSeqObject)`で0.5（50% quantile）が設定されているみたい。`p=0.95`
- quantileにはゼロを除いた値のみが考慮される（quantile(x[x>0], 0.5)）。
- `norm = FALSE`, `log = FALSE`の時は元の値が返ってくる
- `norm = TRUE`での計算は、value*1000/sum(value[value <= quantile])（quantile以下の値の和で割って1000を掛ける）。
- `log = TRUE`での計算は、log2(value + 1)（最低値がゼロになる）もちろんnormの後。

r/CSSnorm.R
```
suppressMessages(library(tidyverse))
suppressMessages(library(metagenomeSeq))


CSSnorm <- function(df){

    df <- as.data.frame(df)

    if(!is.numeric(pull(df,1))){
        col2rowname <- TRUE
        label_colname <- colnames(df)[1]
        rownames(df) <- df[,1]
    } else {
        col2rowname <- FALSE
    }

    df <- select(df, where(is.numeric))

    # convert table into package format
    metaSeqObject <- newMRexperiment(df)
    rm(df); gc(); gc()

    # CSS normalization
    metaSeqObject_CSS <- cumNorm(
      metaSeqObject,
      p = cumNormStatFast(metaSeqObject)
    )
    rm(metaSeqObject); gc(); gc()

    # convert CSS normalized data into data.frame-formatted OTU table (log2 transformed data)
    df_CSS <- data.frame(MRcounts(metaSeqObject_CSS, norm=TRUE, log=TRUE))
    if(isTRUE(col2rowname)){
    df_CSS <- rownames_to_column(df_CSS, label_colname)
    }

    df_CSS <- as_tibble(df_CSS)

    rm(metaSeqObject_CSS); gc(); gc()

    return(df_CSS)

}

```