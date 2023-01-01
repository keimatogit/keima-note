# MLD

[Mean log deviation](https://en.wikipedia.org/wiki/Mean_log_deviation)

```r
mld <- function(value){
    res <- log(mean(value)) - mean(log(value))
    return(res)
}
```

```r
mld_groups <- function(value, group){
    unique_groups <- unique(group)
    share <- map_dbl(unique_groups, sum(group == .x)
    tmp <- map(unique_groups){}
    between <- log(mean(value)) - 
    shares <- map_dbl(unique_groups, ~sum(groups == .x)/length(groups))
    group_means <- map_dbl(unique_groups, ~mean(values[groups == .x]))

    res <- log(mean(values)) - sum(shares*log(group_means))
    return(res)
}
```