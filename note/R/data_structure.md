# オブジェクトのデータ型と内部構造の確認

<!-- @import "[TOC]" {cmd="toc" depthFrom=2 depthTo=6 orderedList=false} -->
<!-- code_chunk_output -->

- [サンプルデータ](#サンプルデータ)
- [データ型の確認](#データ型の確認)
  - [typeof](#typeof)
  - [mode](#mode)
  - [class](#class)
- [内部構造の確認](#内部構造の確認)
  - [str](#str)

<!-- /code_chunk_output -->


## サンプルデータ

```r

# 文字列ベクトル
vec_chr <- c("A", "B", "C")

# 数値ベクトル
vec_int <- 1:3

# 理論値ベクトル
vec_bool <- c(TRUE, FALSE, FALSE)

# 名前つき数値ベクトル
vec_named_int <- 1:3 %>%
  set_names(c("A", "B", "C"))

# 順序付きfactorベクトル
vec_factor <- factor(
  c("A", "B", "C"),
  levels=c("A", "B", "C"),
  ordered = TRUE)

# 行列
mat <- matrix(1:9, nrow=3, ncol=3)

# リスト
lt <- list(
  "obj1" = c("A", "B", "C"),
  "obj2" = 1:2,
  "obj3" = TRUE
  )

# データフレーム
df <- data.frame(
  "col1" = c("A", "B", "C"),
  "col2" = 1:3,
  "col3" = c(TRUE, FALSE, FALSE)
  )

# tibble
tb <- tibble(
  "col1" = c("A", "B", "C"),
  "col2" = 1:3,
  "col3" = c(TRUE, FALSE, FALSE)
  )   

# S4オブジェクト
setClass (
  "person",
  representation (
    name = "character",
    age = "numeric"
  ),
  prototype = list(
    name = NULL,
    age = NULL
  )
)
s4_obj <- new("person", name="keima", age=30)

# 関数
func <- function(x){x*100}
```

## データ型の確認

Rのオブジェクトには、type, mode, classの３種類の型が設定されているらしい。それぞれ`typeof()` `mode()` `class()`で調べる。


### typeof

ベクトルと行列はcharacter, integerなどのデータタイプ（factorはinteger）、リストとデータフレームとtibbleは"list"、S4オブジェクトは"S4"が返ってくる。

```r
# 文字列ベクトル
> typeof(vec_chr)
[1] "character"

# 数値ベクトル
> typeof(vec_int)
[1] "integer"

# 理論値ベクトル
> typeof(vec_bool)
[1] "logical"

# 名前付き数値ベクトル
> typeof(vec_named_int)
[1] "integer"

# 順序付きfactorベクトル
> typeof(vec_factor)
[1] "integer"

# 行列
> typeof(mat)
[1] "integer"

# リスト
> typeof(lt)
[1] "list"

# データフレーム
> typeof(df)
[1] "list"

＃tibble
> typeof(tb)
[1] "list"

# S4オブジェクト
> typeof(s4_obj)
[1] "S4"

# 関数
> typeof(func)
[1] "closure"
```

### mode

ベクトルと行列はcharacter, numericなどのデータタイプ（factorはnumeric）、リストとデータフレームとtibbleは"list"、S4オブジェクトは"S4"が返ってくる。typeとmodeの違いがよく分からない。

```r
# 文字列ベクトル
> mode(vec_chr)
[1] "character"

# 数値ベクトル
> mode(vec_int)
[1] "numeric"

# 理論値ベクトル
> mode(vec_bool)
[1] "logical"

# 名前つき数値ベクトル
> mode(vec_named_int)
[1] "numeric"

# 順序付きfactorベクトル
> mode(vec_factor)
[1] "numeric"

# 行列
> mode(mat)
[1] "numeric"

# リスト
> mode(lt)
[1] "list"

# データフレーム
> mode(df)
[1] "list"

# tibble
> mode(tb)
[1] "list"

# S4オブジェクト
> mode(s4_obj)
[1] "S4"

# 関数
> mode(func)
[1] "function"
```

### class

ベクトルはcharacter、integerなどのデータタイプが返ってくる（factorは"factor"）。行列、リスト、データフレーム、tibbleは、"matrix"、"list"、"data.frame"、"tbl_df"などが返ってくる。S4オブジェクトは定義したクラス名が返ってくる。

```r
# 文字列ベクトル
> class(vec_chr)
[1] "character"

# 数値ベクトル
> class(vec_int)
[1] "integer"

# 理論値ベクトル
> class(vec_bool)
[1] "logical"

# 名前付き数値ベクトル
> class(vec_named_int)
[1] "integer"

# 順序付きfactorベクトル
> class(vec_factor)
[1] "ordered" "factor"

# 行列
> class(mat)
[1] "matrix"

# リスト
> class(lt)
[1] "list"

# データフレーム
> class(df)
[1] "data.frame"

# tibble
> class(tb)
[1] "tbl_df"     "tbl"        "data.frame"

# S4オブジェクト
> class(s4_obj)
[1] "person"
attr(,"package")
[1] ".GlobalEnv"

# 関数
> class(func)
[1] "function"
```

## 内部構造の確認

### str

すごく便利。

```r
# 文字列ベクトル
> str(vec_chr)
 chr [1:3] "A" "B" "C"

# 数値ベクトル
> str(vec_int)
 int [1:3] 1 2 3

# 理論値ベクトル
> str(vec_bool)
 logi [1:3] TRUE FALSE FALSE

# 名前付き数値ベクトル
> str(vec_named_int)
 Named int [1:3] 1 2 3
 - attr(*, "names")= chr [1:3] "A" "B" "C"

# 順序付きfactorベクトル
> str(vec_factor)
 Ord.factor w/ 3 levels "A"<"B"<"C": 1 2 3

# リスト
> str(lt)
List of 3
 $ obj1: chr [1:3] "A" "B" "C"
 $ obj2: int [1:2] 1 2
 $ obj3: logi TRUE

# データフレーム
> str(df)
'data.frame':	3 obs. of  3 variables:
 $ col1: Factor w/ 3 levels "A","B","C": 1 2 3
 $ col2: int  1 2 3
 $ col3: logi  TRUE FALSE FALSE

# tibble
> str(tb)
Classes ‘tbl_df’, ‘tbl’ and 'data.frame':	3 obs. of  3 variables:
 $ col1: chr  "A" "B" "C"
 $ col2: int  1 2 3
 $ col3: logi  TRUE FALSE FALSE

# マトリックス
> str(mat)
 int [1:3, 1:3] 1 2 3 4 5 6 7 8 9

# S4オブジェクト
> str(s4_obj)
Formal class 'person' [package ".GlobalEnv"] with 2 slots
  ..@ name: chr "keima"
  ..@ age : num 30

# 関数
> str(func)
function (x)  
 - attr(*, "srcref")= 'srcref' int [1:8] 1 9 1 26 9 26 1 1
  ..- attr(*, "srcfile")=Classes 'srcfilecopy', 'srcfile' <environment: 0x7fa605d54d98>
```
