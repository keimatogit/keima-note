# Alpha diversity

[α多様性とβ多様性](https://qiita.com/keisuke-ota/items/9701c583e153df467c05)


## Shannon Index

- (割合) x log(割合) の和 x - 1。割合の対数で重みづけされるので、希少な種の重みがちょっと増す。



## Simpson Index

1 - (割合の二乗和)。それ自体の頻度で重みづけされ希少な種に重みるので、主要な種に重みづけされる（希少な種の影響を受けにくい）。多くの種が同じような割合で存在するときに１に近づく。

```r
simpson <- 1 - sum(x^2)
```

vegan::diversityを使う場合
```r
vegan::diversity(x, index = "simpson")
```