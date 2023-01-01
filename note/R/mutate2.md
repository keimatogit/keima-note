# mutateでも使える条件分岐（if_else, case_when, recode, factor）


<!-- @import "[TOC]" {cmd="toc" depthFrom=2 depthTo=6 orderedList=false} -->
<!-- code_chunk_output -->

- [if_else](#if_else)
- [case_when](#case_when)
- [recode/recode_factor](#recoderecode_factor)
- [おまけ：factor関数での順序付き尺度のラベルづけ](#おまけfactor関数での順序付き尺度のラベルづけ)

<!-- /code_chunk_output -->


## if_else

[Vectorised if](https://dplyr.tidyverse.org/reference/if_else.html)

Usage
```r
dplyr::if_else(condition, true, false, missing = NULL)
```

Examples
```r
tibble(sleep) %>%
  mutate(extra.2 = if_else(extra >= 0, "plus", "minus"))

  # # A tibble: 20 x 4
  #    extra group ID    extra.2
  #    <dbl> <fct> <fct> <chr>  
  #  1   0.7 1     1     plus   
  #  2  -1.6 1     2     minus  
  #  3  -0.2 1     3     minus  
  #  4  -1.2 1     4     minus  
  #  5  -0.1 1     5     minus  
  #  6   3.4 1     6     plus   
  #  7   3.7 1     7     plus   
  #  8   0.8 1     8     plus   
  #  9   0   1     9     plus   
  # 10   2   1     10    plus   
  # 11   1.9 2     1     plus   
  # 12   0.8 2     2     plus   
  # 13   1.1 2     3     plus   
  # 14   0.1 2     4     plus   
  # 15  -0.1 2     5     minus  
  # 16   4.4 2     6     plus   
  # 17   5.5 2     7     plus   
  # 18   1.6 2     8     plus   
  # 19   4.6 2     9     plus   
  # 20   3.4 2     10    plus   
```

mutateじゃなくても使える
```r
x <- c(-5:5)
if_else(x >= 0, "plus", "minus")

# [1] "minus" "minus" "minus" "minus" "minus" "plus"  "plus"  "plus"  "plus"
# [10] "plus"  "plus"
```


## case_when

[A general vectorised if](https://dplyr.tidyverse.org/reference/case_when.html)

Usage
```r
dplyr::case_when(...)
```

条件式を順番にみていって、あてはまった時点の値を入れる。最後の`TRUE ~ `は、それまでのどの条件にも当てはまらなかった場合に入れる値。

Examples
```r
tibble(sleep) %>%
  mutate(extra.2 = case_when(
    extra > -2 & extra <= 0 ~ "(-2, 0]",
    extra <= 2              ~ "(0,  2]",
    extra <= 4              ~ "(2,  4]",
    extra <= 6              ~ "(4,  6]",
    TRUE ~ "out of range"
    )
  )

  # # A tibble: 20 x 4
  #    extra group ID    extra.2
  #    <dbl> <fct> <fct> <chr>  
  #  1   0.7 1     1     (0,  2]
  #  2  -1.6 1     2     (-2, 0]
  #  3  -0.2 1     3     (-2, 0]
  #  4  -1.2 1     4     (-2, 0]
  #  5  -0.1 1     5     (-2, 0]
  #  6   3.4 1     6     (0,  4]
  #  7   3.7 1     7     (0,  4]
  #  8   0.8 1     8     (0,  2]
  #  9   0   1     9     (-2, 0]
  # 10   2   1     10    (0,  2]
  # 11   1.9 2     1     (0,  2]
  # 12   0.8 2     2     (0,  2]
  # 13   1.1 2     3     (0,  2]
  # 14   0.1 2     4     (0,  2]
  # 15  -0.1 2     5     (-2, 0]
  # 16   4.4 2     6     (0,  6]
  # 17   5.5 2     7     (0,  6]
  # 18   1.6 2     8     (0,  2]
  # 19   4.6 2     9     (0,  6]
  # 20   3.4 2     10    (0,  4]
```

```r
x <- c(1:10)
case_when(
  x %% 2 == 0 ~ "even",
  x %% 3 == 0 ~ "multiple of 3",
  x %% 5 == 0 ~ "multiple of 5",
  TRUE ~ "others"
)

# [1] "others"        "even"          "multiple of 3" "even"         
# [5] "multiple of 5" "even"          "others"        "even"         
# [9] "multiple of 3" "even"             
```

## recode/recode_factor

[Recode values](https://dplyr.tidyverse.org/reference/recode.html)

Usage

```r
recode(.x, ..., .default = NULL, .missing = NULL)
recode_factor(.x, ..., .default = NULL, .missing = NULL, .ordered = FALSE)
```
Examples
```r
tibble(sleep) %>%
  mutate(group.2 = recode(group, `1` = "drug A", `2` = "drug B"))

# # A tibble: 20 x 4
#    extra group ID    group.2
#    <dbl> <fct> <fct> <fct>  
#  1   0.7 1     1     drug A
#  2  -1.6 1     2     drug A
#  3  -0.2 1     3     drug A
#  4  -1.2 1     4     drug A
#  5  -0.1 1     5     drug A
#  6   3.4 1     6     drug A
#  7   3.7 1     7     drug A
#  8   0.8 1     8     drug A
#  9   0   1     9     drug A
# 10   2   1     10    drug A
# 11   1.9 2     1     drug B
# 12   0.8 2     2     drug B
# 13   1.1 2     3     drug B
# 14   0.1 2     4     drug B
# 15  -0.1 2     5     drug B
# 16   4.4 2     6     drug B
# 17   5.5 2     7     drug B
# 18   1.6 2     8     drug B
# 19   4.6 2     9     drug B
# 20   3.4 2     10    drug B
```

named vectorを使うこともできる
```r
replacement <- c("drug A", "drug B") %>%
  set_names(c(1, 2))

tibble(sleep) %>%
  mutate(group.2 = recode(group, !!!replacement))

# # A tibble: 20 x 4
#    extra group ID    group.2
#    <dbl> <fct> <fct> <fct>  
#  1   0.7 1     1     drug A
#  2  -1.6 1     2     drug A
#  3  -0.2 1     3     drug A
#  4  -1.2 1     4     drug A
#  5  -0.1 1     5     drug A
#  6   3.4 1     6     drug A
#  7   3.7 1     7     drug A
#  8   0.8 1     8     drug A
#  9   0   1     9     drug A
# 10   2   1     10    drug A
# 11   1.9 2     1     drug B
# 12   0.8 2     2     drug B
# 13   1.1 2     3     drug B
# 14   0.1 2     4     drug B
# 15  -0.1 2     5     drug B
# 16   4.4 2     6     drug B
# 17   5.5 2     7     drug B
# 18   1.6 2     8     drug B
# 19   4.6 2     9     drug B
# 20   3.4 2     10    drug B
```

recodeをfactorに使用すると、元のlevelsが保たれる。
recode_factorを使用すると、levelsが更新されreplacementの順番になる。

```r
tibble(sleep) %>%
  mutate(group.2 = recode(group, `2` = "drug B", `1` = "drug A")) %>%
  pull(group.2) %>%
  levels
# [1] "drug A" "drug B"

tibble(sleep) %>%
  mutate(group.2 = recode_factor(group, `2` = "drug B", `1` = "drug A")) %>%
  pull(group.2) %>%
  levels
# [1] "drug B" "drug A"
```

## おまけ：factor関数での順序付き尺度のラベルづけ

levelsに対してlabelsを貼れる。

```r
levels(tibble(sleep)$group)
# [1] "1" "2"

tibble(sleep) %>%
  mutate(group.2 = factor(group, labels=c("drug A", "drug B")))
# # A tibble: 20 x 4
#    extra group ID    group.2
#    <dbl> <fct> <fct> <fct>  
#  1   0.7 1     1     drug A
#  2  -1.6 1     2     drug A
#  3  -0.2 1     3     drug A
#  4  -1.2 1     4     drug A
#  5  -0.1 1     5     drug A
#  6   3.4 1     6     drug A
#  7   3.7 1     7     drug A
#  8   0.8 1     8     drug A
#  9   0   1     9     drug A
# 10   2   1     10    drug A
# 11   1.9 2     1     drug B
# 12   0.8 2     2     drug B
# 13   1.1 2     3     drug B
# 14   0.1 2     4     drug B
# 15  -0.1 2     5     drug B
# 16   4.4 2     6     drug B
# 17   5.5 2     7     drug B
# 18   1.6 2     8     drug B
# 19   4.6 2     9     drug B
# 20   3.4 2     10    drug B
```

```r
x <- c(1,2,1,1,2,3,3)
factor(x, levels = c(1,2,3), labels = c("low", "mid", "high"))
# [1] low  mid  low  low  mid  high high
# Levels: low mid high
```
