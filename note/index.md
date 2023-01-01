
<!-- @import "[TOC]" {cmd="toc" depthFrom=2 depthTo=3 orderedList=false} -->
<!-- code_chunk_output -->

- [R](#r)
- [python](#python)
- [データベース・SQL](#データベースsql)
- [github](#github)
- [bash](#bash)
- [バイオインフォ](#バイオインフォ)
- [その他](#その他)

<!-- /code_chunk_output -->

## R
- [オブジェクトのデータ型と内部構造の確認(typeof, mode, class, str)](R/data_structure.md)
- [mutateでも使える条件分岐（if_else, case_when, recode, (factor)）](R/mutate2.md)
- [chop and unchop]()
- [nest and unnest]()
- [factor型について](R/factor.md)
- [パッケージの場所（.libPaths()）](R/libPaths.md)
- [renvを使用したパッケージ管理](R/renv.md)
- [パッケージの作り方](R/make_package.md)
- [Rmarkdownの使い方](R/rmarkdown.md)
- [データベースを使用する（RMySQL/MariaDB）](R/database.md)


- [統計量計算用関数いろいろ](R/stat_values.md)

## python
- [Jupyter notebookのいろいろ](python/jupyter.md)
- [データベースを使用する(MariaDB)](python/database.md)
- [seabornサンプルデータの使い方](python/sample_data.md)
- [グラフ](python/graphs.md)
- [ipywidgetsでインタラクティブなグラフを作る](python/ipywidgets.md)

## データベース・SQL

準備とDB操作
- [環境設定（MariaDBのインストールと起動）](sql/mariadb_install.md)
- [DBサーバーへの接続とSQLの実行方法](sql/execute_sql.md)
- [DBユーザーの管理](sql/user.md)
- [データベースの操作](sql/database.md)

テーブルの管理と構造
- [テーブルの操作](sql/table.md)
- [インデックス（キー）](sql/sql_index.md)
- [外部キー制約、check制約](sql/constraint.md)

データの書き換え
- [レコード（行）の追加・編集・削除](sql/record.md)
- [csvファイルのインポート(LOAD DATA INFILE)](sql/csv_import.md)

データの取得
- [データの取得（SELECT）① - 基本 -](sql/select1.md)
- [データの取得（SELECT）② - サブクエリの利用 -](sql/select2.md)
- [テーブルの結合（JOIN）](sql/join.md)

その他
- [MySQLサンプルデータの使い方](sql/sample_data.md)
- [SQLiteの使い方](sql/sqlite.md)

## github
- [githubの使い方](github.md)
- [Guthub Pagesの使い方](github_pages.md)


## bash

## バイオインフォ
- FastQC
- seqkit
- MEGAHIT
- MetagenMark
- CD-HIT
- Bowtie2
- samtools
- coverM
- HUMAnN3

## その他
- [コマンドラインでのプログラム実行時に引数を扱う（python, R, shell）](commandargs.md)
- [condaの使い方](conda.md)
- [singularityの使い方](github.md)
- [WindowsでUbuntuをを使う](ubuntu_setup.md)
- [ssh接続でホスト名を登録](ssh_hostname.md)
- [テキストエディタATOM](atom.md)
