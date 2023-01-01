# renvを使用したパッケージ管理

参考サイト：
[Introduction to renv](https://rstudio.github.io/renv/articles/renv.html)  
[Collaborating with renv](https://rstudio.github.io/renv/articles/collaborating.html)  
[Rのパッケージ管理のためのrenvの使い方](https://qiita.com/okiyuki99/items/688a00ca9a58e42e3bfa)  
[renv: Project Environments for R](https://blog.rstudio.com/2019/11/06/renv-project-environments-for-r/)

```
install.packages("renv")

# バージョン確認
library(renv)
packageVersion("renv") 
```

## project libraryの作り方


1.  `renv::init()` でrenvによるパッケージ管理を始める。Rプロジェクトの場所に、renv/というprivate R libraryの置き場と、renv.lock、.Rprofileが作られる。.Rprofileにはsource("renv/activate.R") が記載されており、Rプロジェクト開始とともにrenvによるプロジェクト管理をスタートする処理が書かれている。

2. Rプロジェクトで必要なパッケージをインストールしたりアンインストールしたり、通常の処理をする

3.  `renv::snapshot()` でproject libraryの状態をlockfile (renv.lock)に保存。インストールしたパッケージ以外にも、たとえばプロジェクト中のスクリプトにlibrary(ggplot2)と書かれていると、ggplot2が使われていることをrenvが発見してrenv.lockに書き込んでくれる。

4. 処理を続け、適宜 `renv::snapshot()` を実行することで、project libraryの状態を保存してくれる。 `renv::restore()` では以前の保存状態に戻すことができる。

5. renvによるパッケージ管理を解除するには `renv::deactivate()`

補足メモ：
- パッケージ情報はrenv.lockに書き込まれ、基本的にこのファイルだけあれば再現できるみたい。`renv::snapshot()`実行後、テキストエディタで開いて書き込まれているか確認してみよう。
-  `renv::init()` でパッケージ管理を始めるとパッケージが何もない状態から始まるので、インストールにけっこう時間がかかる。
-  `renv::init()` は一番初めだけで、あとは `renv::snapshot()` をしていく感じ。
- googleドライブでは自動で作成される"Icon?"ファイル（renv/library/R-4.0/x86_64-apple-darwin17.0/Icon）のせいでそんなパッケージないとmessageが出た。消しとくと良い。


## project libraryの再現

renv.lockと自動読み込みのためのファイル（.Rprofile と renv/activate.R。renvフォルダは丸ごと持っていってしまうと良い）をプロジェクト場所（作業ディレクトリ）にコピーして `renv::restore()`

または

renv.lockのみコピーして `renv::init()`  `renv::restore()`を実行（他のはinit()で生成される）。


## 

```
renv::settings$external.libraries(Sys.getenv("R_LIBS_USER"))
```