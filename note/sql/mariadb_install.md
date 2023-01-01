# 環境設定（MariaDBのインストールと起動）

リレーショナルデータベース管理システム(RDBMS)には、MySQL、MariaDB、SQliteなど色々な製品がありますが、ここではMySQL派生のMariaDBを使用しています。
データベースを操作するSQL言語は基本的にどのDBMSでも共通ですが、DBMSによってそのシステム専用のSQLも存在します。


<!-- @import "[TOC]" {cmd="toc" depthFrom=2 depthTo=6 orderedList=false} -->
<!-- code_chunk_output -->

- [MariaDBのインストール(Mac)](#mariadbのインストールmac)
- [MariaDBの起動](#mariadbの起動)
- [パケット通信制限の確認と変更](#パケット通信制限の確認と変更)

<!-- /code_chunk_output -->


## MariaDBのインストール(Mac)

Homebrewでインストール。
```
brew update                 # Homebrewのアップデート
brew install mariadb        # MariaDBのインストール
sudo chown -R user_name  /usr/local/Homebrew #パーミッション断られたとき
brew services start mariadb # MariaDBの起動
mysql_secure_installation   # 初期設定
```

## MariaDBの起動

brewでインストールした場合
```
brew services start mariadb   # 起動
brew services restart mariadb # 再起動
brew services stop mariadb    # 終了
status                        # ステータス確認
```

通常はこう？
```
mysql.server start   # 起動
mysql.server restart # 再起動
mysql.server stop    # 終了
mysql.server status  # ステータス確認
```

## パケット通信制限の確認と変更

[max_allowed_packetの設定変更](https://hodalog.com/how-to-set-max-allowed-packet/)

通信制限量の確認（デフォルトは16MB？）
```
SHOW VARIABLES LIKE 'max_allowed_packet';
```
変更するには、設定ファイルに`max_allowed_packet = 64M`等を追加する。（`mysql --help`で設定ファイルの場所が分かるので、そこでインクルードしてるファイルを見て該当箇所を書き換える。`/etc/mysql/conf.d/mysqldump.cnf`など？）それかユーザー用設定フォルダにファイルを差s区政して書き込む`~/.my.cnf/my_mysql_setting.cnf`など。
