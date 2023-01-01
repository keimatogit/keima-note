# パッケージの場所

[libPaths {base}](https://stat.ethz.ch/R-manual/R-devel/library/base/html/libPaths.html)
[libPaths: Search Paths for Packages](https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/libPaths)

## パッケージの場所一覧

Rが認識しているパッケージの場所一覧を表示します。

```
.libPaths()
```

## インストール時

デフォルトでは`.libPaths()`で初めに出力される場所に保存されます。他の場所に保存する場合は`lib = "" `で指定します。
```
install.packages("tidyverse") # .libPaths()で初めに表示される場所に保存される
install.packages("tidyverse", lib = "./library") # 場所を指定して保存
```

## 読み込み時

デフォルトでは`.libPaths()`で表示される場所を探して読み込みます（.libPath()の出力順に探すので、同じパッケージが複数のフォルダにある場合は先に表示されるフォルダのものが読み込まれます）。他の場所にあるパッケージを読み込む場合は`lib.loc = "" `で指定します。

```
library(tidyverse)
library(tidyverse, lib.loc = "./library")
```

## .libPathの変更

パッケージのインストールや読み込み場所を変更したい場合、`libPath()`に指定したい場所のパスを渡すと、.libPath()の中身が、指定パス、R_LIBS_SITE（ない場合もある）、システムのパスのみに変更されます。

```
.libPath("./library")
```
複数指定したい場合はパスをベクトルで渡します。追加したい場合は.libPath()を含めたベクトルを渡します。

```
.libPath(c(.libPaths(), "./library"))
```

## いつも変える場合

.Renvironで環境変数を設定しておくといい。
- R_LIBS: プロジェクトごとの一時的な設定？読み書き時にR_LIBS_USERよりも優先される。
- R_LIBS_USER: ユーザーが使いたいライブラリ。他のユーザーには影響しない。R_LIBS_SITEよりも優先される。
- R_LIBS_SITE: あるサイト（例えば、同じコンピュータ）で複数のユーザーが共有するライブラリ。書き込みはできない？

以下の変換指定子を使用できる。
```
%V
R version number including the patchlevel (e.g., 2.5.0).

%v
R version number excluding the patchlevel (e.g., 2.5).

%p
the platform for which R was built, the value of R.version$platform.

%o
the underlying operating system, the value of R.version$os.

%a
the architecture (CPU) R was built on/for, the value of R.version$arch.
```