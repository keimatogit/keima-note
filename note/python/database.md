# データベースを使用する

- pandas.DataFrameの`read_sql`と `to_sql`が使いやすい。sqlalchemy + mysql-connector-pythonでの接続はどちらも使用できたが、mariadbでの接続ではto_sqlが実行できなかったので注意。
- データフレームのインポートでは、mariadbのexacutemanyは速いけど使いにくい。to_sqlは使いやすいけどやや遅い（RのRMariaDBのdbWriteTableより遅い）。どうしても速くしたいのでなければto_sqlが良い。

----
TOC
<!-- @import "[TOC]" {cmd="toc" depthFrom=2 depthTo=6 orderedList=false} -->
<!-- code_chunk_output -->

- [sqlalchemy + mysql-connector-pythonで接続](#sqlalchemy-mysql-connector-pythonで接続)
- [mariadbライブラリで接続](#mariadbライブラリで接続)
- [DataFrame.to_sqlの挙動について](#dataframeto_sqlの挙動について)

<!-- /code_chunk_output -->


## sqlalchemy + mysql-connector-pythonで接続
sqlalchemy、mysql-connector-pythonライブラリが必要

```python

from sqlalchemy import create_engine
import pandas as pd

engine = create_engine('mysql+mysqlconnector://user_name:password@host/database_name')

# read_sql(query結果をデータフレームで取得)
query = 'SELECT * FROM table_name'
df = pd.read_sql(query, engine)

# to_sql(データフレームをDBにインポート)
# if_exists{‘fail’, ‘replace’, ‘append’}, default ‘fail’
df.to_sql('table_name', engine, index=False, if_exists='append')
```

## mariadbライブラリで接続
[How to connect Python programs to MariaDB](https://mariadb.com/ja/resources/blog/how-to-connect-python-programs-to-mariadb/)

まずはDBサーバーに接続
```python

import mariadb
import sys
import pandas as pd

# Connect to MariaDB Platform
try:
  conn = mariadb.connect(
    user       = 'user_name',
    password   = 'password',
    dbname     = 'db_name',
    host       = 'host_name',
    autocommit = False # デフォはFalse
    )
except mariadb.Error as e:
    print(f"Error connecting to MariaDB Platform: {e}")
    sys.exit(1)
```

データの取得はpandas.read_sqlを使うと便利
（mariadb.connectの接続ではto_sqlはできなかった。変わるかもしれないのでまた試してみてほしい。）
```python
query = 'SELECT * FROM table_name'
df = pd.read_sql(query, conn)
```

cursorを作ってexecuteでSQLを実行
```sql
# Get Cursor
cur = conn.cursor()

# テーブル作成
cur.execute('CREATE TABLE table_name')
conn.commit()

# レコードの追加（１件）
cur.execute(
    "INSERT INTO table_name (col1, col2, col3) VALUES (?, ?, ?)",
    (value1, value2, value3))
conn.commit()
```

executemanyでのデータフレームのインポート。
（pandasのto_sqlより速かったが、nanがnull認識されないのでNoneにしておく必要がある、'?'のリストが必要、など使いにくい感じなので、よほど急ぎでなければto_sqlの方が良いと思う。）

```python

insert_field_names = ['col1', 'col2', 'col3']
questions = ['?'] * len(insert_field_names)

df = df.where(pd.notnull(df), None) # NaNをNoneに変換

try:
    cur = conn.cursor()
    cur.executemany(
        f'INSERT INTO table_name ({", ".join(insert_field_names)}) VALUES ({", ".join(questions)})',
        [tuple(x) for x in df.values] # データはタプルのリストでないといけない
    )
except mariadb.Error as e:
    print(f'Error: {e}')

conn.commit()
```

最後に接続解除する
```
conn.close()
```

## DataFrame.to_sqlの挙動について

- データフレームのカラム名に、DB側テーブルのフィールド名にないものがある場合、エラーが出て読み込まれなかった。
- DB側テーブルのフィールド名に、データフレームのカラム名にないものがある場合、そのフィールドはNULL（またはデフォルト値）が入力された。
- DBテーブルの制約（check制約、FK制約など）を満たさない行が1行でもあると、データフレーム全てインポートされなかった。
