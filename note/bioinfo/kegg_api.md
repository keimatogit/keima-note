# KEGG API

[KEGG API](https://www.kegg.jp/kegg/docs/keggapi.html)


KEGG KO（やpathwayやmodule）のidから名前（definition）を取得する関数例

```
library(tidyverse)

getKEGGnames <- function(ids, db="ko"){
    url <- paste0("https://rest.kegg.jp/list/", db)
    db <- read_tsv(url, col_types=cols(), col_names = c("id", "name"))
    db <- db %>%
        mutate(id = str_remove(id, "^ko:|^path:|md:")) # ひとまず取り除くことにしておく
    
    names <- db$name[match(ids, db$id)]
    
    return(names)

}

```


KEGG koからpathwayやmoduleを取得する関数。
複数pathway/moduleに当たるので、基本的に重複を除いたデータフレーム（変換元・変換先の２列）で返す。left_joinで使うことを想定。
return_list=TRUEでリスト形式で返す（時間がかかるかも）。

```
library(tidyverse)

linkKEGGids <- function(ids, from="ko", to="pathway", return_list=FALSE){

    url <- paste0("https://rest.kegg.jp/link/", to, "/", from)
    db <- read_tsv(url, col_types=cols(), col_names = c("from", "to"))
    db <- db %>%
        filter(if_all(everything(), ~!str_detect(., "^path:ko"))) # pathwayはpath:koとpath:mapがある。list/pathwayはpath:mapなのでそちらに揃えておく
    db <- db %>%
        mutate(across(everything(), ~str_remove(., "^ko:|^path:|md:"))) # ひとまず取り除くことにしておく
    db <- db %>%
        filter(from %in% unique(ids)) # データを減らしておく
    
    if (isTRUE(return_list)){
        res <- map(ids, ~filter(db, from == .x) %>% pull(to))
        
    } else {
       res <- db %>%
            mutate(from = factor(from, levels=unique(ids))) %>% # 入力順に並べ替えておく
            arrange(from) %>%
            mutate(from = as.character(from)) %>%
            rename(c("from" ,"to") %>% set_names(from, to))
    }
    
    return(res)
    
}
```
