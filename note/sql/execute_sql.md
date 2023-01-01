# DBサーバーへの接続とSQLの実行方法

データベース(DB)サーバー内に複数の「データベース」があり、各「データベース」に複数のテーブルが格納されています。SQL実行の際はDBサーバーだけでなく使用するデータベース名を指定する必要があります。


<!-- @import "[TOC]" {cmd="toc" depthFrom=2 depthTo=6 orderedList=false} -->
<!-- code_chunk_output -->

- [基本：対話モード](#基本対話モード)
- [スクリプトファイルで実行](#スクリプトファイルで実行)
- [ワンライナーで実行](#ワンライナーで実行)
- [python、Rから](#python-rから)

<!-- /code_chunk_output -->


## 基本：対話モード

データベース名を指定してDBサーバーに接続
- -pはパスワード要求。ユーザーにパスワードが設定されている場合必要。
- ユーザー名のデフォルトはPCのログインユーザー名。
- ホスト名のデフォルトはlocalhost。
```
mysql database_name [-u user_name] [-h host_name] [-p] [-P port_number]
```

またはDBサーバーに接続後、使用データベース名を指定
```
mysql [-u user_name] [-h host_name] [-p]
use database_name;
```

以下のように表示されていたらDBに接続している状態なので、ここでSQLを使える。
```
MariaDB [(database_name)]>
```
終るときは、`Ctrl+C`で接続解除。

ちなみに対話モード中にsqlファイルを実行する場合は、`source [スクリプト名]`。
```sql
source script.sql;
```

## スクリプトファイルで実行

SQL文をまとめて書いて保存しておき、それを読み込んでを実行することもできる。例えば以下のようなexample.sqlを作成する。
```sql
# example.sql

SELECT user_name AS name, count(*) AS N from user_table
GROUP BY age
ORDERED BY N desc;
```

DBサーバーに接続する際のコマンドの最後に`< スクリプト名`をつけると、対象DBでスクリプト内のSQLを実行してくれる。
```
mysql database_name [-u user_name] [-h host_name] [-p] < example.sql
```

さらに`> [保存先]`をつけると、実行結果を直接表示せず.txtや.csvファイルとして保存できる。デフォルトタブ区切りなので、カンマ区切りなどにする際はsedで置き換える。
```
mysql database_name [-u user_name] [-h host_name] [-p] < example.sql > result.csv

# カンマ区切りにするとき
mysql database_name [-u user_name] [-h host_name] [-p] < example.sql | sed -e 's/\t/,/g' > result.txt
```

`<`で読み込むほかにも、`cat`や`-e`オプションを使う方法もある。
```
cat example.sql | mysql DATABASE
mysql testdb -e "source example.sql;"
```

## ワンライナーで実行

簡単な指示を出す場合などは、DBサーバーに接続する際のコマンドの最後に`-e "SQL文;"`をつけると、対象DBでそのSQLを実行してくれる。
```
mysql database_name [-u user_name] [-h host_name] [-p] -e "show tables;"
```

さらに`> [保存先]`をつけると、実行結果を直接表示せず.txtや.csvファイルとして保存できる。
```
mysql database_name [-u user_name] [-h host_name] [-p] -e "show tables;" > result.txt

# カンマ区切りにするとき
mysql database_name [-u user_name] [-h host_name] [-p] -e "show tables;" | sed -e 's/\t/,/g' > result.txt
```

`-e`オプションを使うほかにも、`echo`や`<<<`を使う方法もある。
```
echo "show tables;" | mysql database_name
mysql database_name <<<'show tables;'
```
## python、Rから


 python: [データベースを使用する(MariaDB)](../python/database.md)
 R: [データベースを使用する(MariaDB)](../R/database.md)
