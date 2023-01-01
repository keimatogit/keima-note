# csvを標準出力に出す

```
if(!is.na(args$modified_clstr)){
        if(args$modified_clstr == "stdout"){
            cat(format_tsv(res$modified_clstr)) # ここ
        } else {
        write_tsv(res$modified_clstr, args$modified_clstr)
        }
    }
```