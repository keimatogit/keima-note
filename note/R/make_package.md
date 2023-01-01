# パッケージの作り方

参考サイト：
[WRITING R PACKAGES IN RSTUDIO](https://ourcodingclub.github.io/tutorials/writing-r-package/) こちらに沿って行います
[devtools — Rパッケージ作成支援](https://heavywatal.github.io/rstats/devtools.html)


<!-- @import "[TOC]" {cmd="toc" depthFrom=2 depthTo=6 orderedList=false} -->
<!-- code_chunk_output -->

- [パッケージ作成](#パッケージ作成)
  - [準備](#準備)
  - [Rスクリプトを用意](#rスクリプトを用意)
  - [DESCRIPTION作成](#description作成)
  - [ヘルプメッセージを用意](#ヘルプメッセージを用意)
  - [データを用意](#データを用意)
- [pkgdownでwebサイトを作る](#pkgdownでwebサイトを作る)
- [.Rbuildignoreについて](#rbuildignoreについて)
- [githubに上げてみよう](#githubに上げてみよう)

<!-- /code_chunk_output -->

## パッケージ作成

### 準備

必要パッケージ
```
install.packages("devtools")
install.packages("roxygen2")
```

パッケージ全体を格納するフォルダ（ここではpackage_dir）を作成し、この中で作業します。package_dir直下にRStudioのFile > New Project...でpackage_dirにR projectを作っておくと便利です。
```
mkdir package_dir
```

### Rスクリプトを用意

package_dirの中に「R」フォルダを作り、ここに関数を書いた「.R」スクリプトファイルを入れます。（「R」フォルダにはいくつでも「.R」ファイルを入れることができるし、各ファイルにはいくつでも関数を入れることができる。）

たとえばこちらの「temp_conversion.R」を「R」フォルダに入れます。
```r
F_to_C <- function(F_temp){
    C_temp <- (F_temp - 32) * 5/9;
    return(C_temp);
}

C_to_F <- function(C_temp){
    F_temp <- (C_temp * 9/5) + 32;
    return(F_temp);
}
```
### DESCRIPTION作成

package_dir直下(Rフォルダの外側)に拡張子なしのDESCRIPTIONというプレーンテキストファイルを作成し、Rパッケージのメタデータを記述します。 `usethis::use_description()`で生成することもできます。
メタデータの書き方：[Package metadata](https://r-pkgs.org/description.html)

```
# シンプルなdescriptionの例
Package: SCCTempConverter
Title: Temperature Conversion Package for Demonstration
Version: 0.0.1.0
Description: What the package does (one paragraph)
Encoding: UTF-8
```
これでパッケージの概形ができました。package_dir直下で`devtools::load_all(".")`を実行し、パッケージが読み込まれて関数が使用できるか確認しましょう（DESCRIPTIONの記述方法が間違っていると注意されます）。

```
devtools::load_all(".")
F_to_C(79)
C_to_F(20)
```

### ヘルプメッセージを用意

`help(F_to_C)`や`?F_to_C`でヘルプメッセージを表示するには、「.R」ファイル内の関数に情報を追記します。「temp_conversion.R」では下のように書かれています。
書き方：[Introduction to roxygen2](https://cran.r-project.org/web/packages/roxygen2/vignettes/roxygen2.html)

```
#' Fahrenheit conversion
#'
#' Convert degrees Fahrenheit temperatures to degrees Celsius
#' @param F_temp The temperature in degrees Fahrenheit
#' @return The temperature in degrees Celsius
#' @examples
#' temp1 <- F_to_C(50);
#' temp2 <- F_to_C( c(50, 63, 23) );
#' @export
F_to_C <- function(F_temp){
    C_temp <- (F_temp - 32) * 5/9;
    return(C_temp);
}

#' Celsius conversion
#'
#' Convert degrees Celsius temperatures to degrees Fahrenheit
#' @param C_temp The temperature in degrees Celsius
#' @return The temperature in degrees Fahrenheit
#' @examples
#' temp1 <- C_to_F(22);
#' temp2 <- C_to_F( c(-2, 12, 23) );
#' @export
C_to_F <- function(C_temp){
    F_temp <- (C_temp * 9/5) + 32;
    return(F_temp);
}
```
project_dir直下で`roxygen2::roxygenise()`を実行すると、これらの情報をもとに「man」フォルダと「NAMESPACE」ファイルが作成されます。

パッケージをロードしてヘルプが表示されることを確認しましょう。
```
devtools::load_all(".")
?F_to_C
```
メッセージを変更したい場合は、.Rファイルを書き換えて `roxygen2::roxygenise()` で更新します（「man」「NAMESPACE」を自分で書き換えない）。


### データを用意

パッケージで呼び出したいデータがある場合、package_dir内に「data」フォルダを作り、ここにデータフレームやベクトルなどのデータを格納した「.rda」ファイルを入れます。

[save: Save R Objects](https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/save)

「.rda」にはいくつでもデータが入りますが、基本的にはひとつのデータをひとつの.rdaファイルにして、データの名前と.rdaファイル名を揃えておくのが分かりやすいかと思います。

```r
Tokyo_2020 <-data.frame(
  month = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
            "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"),
  mean = c(7.1, 8.3, 10.7, 12.8, 19.5, 23.2,
           24.3, 29.1, 24.2, 17.5, 14.0, 7.7)
)

save(Tokyo_2020, file = "data/Tokyo_2020.rda")
```

パッケージをロードしてデータが読み込めるか確認しましょう。`data()`で.rdaファイル名を指定して中のデータを読み込みます。
```
devtools::load_all(".")
data(Tokyo_2020)
```

複数のデータをひとつの.rdaファイルに格納した場合は、全てのデータが読み込まれます。

```r
Tokyo_2019 <-data.frame(
  month = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
            "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"),
  mean = c(5.6, 7.2, 10.6, 13.6, 20.0, 21.8,
           24.1, 28.4, 25.1, 19.4, 13.1, 8.5)
)

Tokyo_2020 <-data.frame(
  month = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
            "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"),
  mean = c(7.1, 8.3, 10.7, 12.8, 19.5, 23.2,
           24.3, 29.1, 24.2, 17.5, 14.0, 7.7)
)

save(Tokyo_2019, Tokyo_2020, file = "data/Tokyo.rda")
```
```
devtools::load_all(".")
data(Tokyo) # Tokyo_2019とTokyo_2020が読み込まれる
```

## pkgdownでwebサイトを作る

[pkgdown](https://pkgdown.r-lib.org/index.html): パッケージ用webサイトを簡単に作ることができます。

```
install.packages("pkgdown")
```

DESCRIPTIONファイルに情報が揃っていないとエラーになるので（著者などもちゃんと書かないといけない）まずはこれをきちんと書きましょう。 `usethis::use_description()`で見本を生成して書き換えると良いかも。

README.mdがホームとなるので、ない場合は適当に作成しておきます。

pkgdownの使用を開始します。「_pkgdown.yml」が作成され、「.Rbuildignore」を書き換えてくれます（なかった場合作成されます）。
```
usethis::use_pkgdown()
```
webサイトを作成します。「docs」フォルダが作成され、中にwebサイトができます。
```
pkgdown::build_site()
```
## .Rbuildignoreについて

「.Rproj」などパッケージの作成に使ったけどパッケージを使うときにいらないものは、project_dir直下の「.Rbuildignore」に書いておきます。以下はpkgdown（後述）が書いてくれた.Rbuildignoreです。

```
^.*\.Rproj$
^\.Rproj\.user$
^_pkgdown\.yml$
^docs$
^pkgdown$
```
## githubに上げてみよう

project_dir内をgithubにアップすると、インストールして使用できます。

githubで公開レポジトリを作成してから、ターミナルでアップロード（[githubの使い方](../github.md)）
```
cd package_dir
git init
git add .
git commit -m "first commit"
git remote add origin repository_url
git push -u origin main
```
Rでインストール＆読み込み

```r
devtools::install_github("github_id/repository_name");
library(SCCTempConverter);
F_to_C(30);
```
