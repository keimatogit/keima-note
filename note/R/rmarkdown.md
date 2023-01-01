# RMarkdonwの使い方

（RstudioではFile -> Nwe file -> Rmarkdownで見本のファイルを作成してくれるので、そこから書き換えるとよいです。）

<!-- @import "[TOC]" {cmd="toc" depthFrom=2 depthTo=6 orderedList=false} -->
<!-- code_chunk_output -->

- [必要なもの](#必要なもの)
- [YAMLヘッダー](#yamlヘッダー)
  - [HTML出力のヘッダー](#html出力のヘッダー)
    - [目次(toc)をつける](#目次tocをつける)
    - [データフレームの表示方法](#データフレームの表示方法)
- [チャンクオプション](#チャンクオプション)
  - [グローバル設定](#グローバル設定)
  - [各コードチャンクの設定](#各コードチャンクの設定)

<!-- /code_chunk_output -->


## 必要なもの

RStudioを使用している場合はすでに入っているはず。

Rパッケージ
- rmarkdown
- knitr

ソフト
- [pandoc](https://pandoc.org/installing.html)

## YAMLヘッダー

### HTML出力のヘッダー

[HTML document](https://bookdown.org/yihui/rmarkdown/html-document.html)
[R Markdownのhtml_documentで指定できるyamlヘッダ項目について](https://qiita.com/kazutan/items/726e03dfcef1615ae999)

```
---
title: "my report"
author: "keima"
date: "2021/06/13"
abstract: "abstract of this document."
output: html_document
---
```

#### 目次(toc)をつける

- `toc: true`：tocをつける
- `toc_depth: [integer]`：どの見出しレベルまでtocをつけるか(3だと # title から ### title まで)
- `number_section: true`：通し番号をつける
- `toc_float: true`：tocをサイドに浮かせて表示させる（デフォルトは浮かせず上部に表示）`collapsed: false`で細かいレベルの見出しを畳まずに表示（デフォルトは畳んで表示）

```
---
title: "my report"
author: "keima"
date: "2021/06/13"
output:
  html_document:
    toc: true
    toc_depth: 2
    number_section: true
    toc_float:
      collapsed: false
---
```

#### データフレームの表示方法

- `df_print: "[option]"`：データフレームの表示形式を指定

```
---
title: "my report"
author: "keima"
date: "2021/06/13"
output:
  html_document:
    df_print: "kable"
---
```

## チャンクオプション

[Code Chunks](https://rmarkdown.rstudio.com/lesson-3.html)

- `echo = TRUE`: コード自体を表示するか
- `include = TRUE`: コード自体とコードによる出力を表示するか
- `eval = TRUE`: 実際にコードとして評価（実行）するか

- `warning = TRUE`: Warningを表示するか
- `error = TRUE`: Errorを表示するか
- `message = TRUE`: message関数の結果を表示するか（パッケージ読み込み時など）
- `warning = TRUE`: プロンプトの文字(例：>)をコード実行結果に追加するか

- `fig.width = 7`: 図の幅
- `fig.height = 7`: 図の高さ
- `fig.cap = ""`: 図のキャプション

### グローバル設定

冒頭で`knitr::opts_chunk$set`で設定。

例）コード部分を表示しない＆メッセージを表示しない場合
`````r
```{r setup, include=FALSE}
knitr::opts_chunk$set(
  echo = FALSE,
  message = FALSE
)
```
`````

### 各コードチャンクの設定

`{r }`内に書く。
`````r
```{r chunk_name, chank_option1, chunk_option2, ...}

```
`````

例）パッケージ読み込みなどの準備のとき
`````r
```{r package, include=FALSE}
library(tidyverse)
```
`````
