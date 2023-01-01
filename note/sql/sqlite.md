# SQLiteの使い方

- SQLiteは、MySQLやMariaDBのようなサーバーではなくライブラリとして動くので、よりコンパクトで気軽に使えます。
- SQLiteのデータベースはひとつのファイルとして保存されます。ファイル名はなんでも大丈夫ですが、分かりやすいように拡張子は`.sqlite`にしておきましょう。
- データベースの読み込みや書き込みの権限は、DBファイル権限がそのまま適応されます。 `chmod`で設定しましょう。


## インストール（conda）

conda環境には既に入っていることが多いと思いますが、なければインストールします。pythonから使用するだけなら、pythonの標準ライブラリsqlite3だけあれば大丈夫です。

```bash
conda install -c anaconda sqlite
which sqlite3 # 確認
```

## 起動と終了

[使用DBファイルを指定して] 開始（指定したファイルがない場合は作成されます）

|   |   |
|---|---|
| [使用DBファイルを指定して]開始  |  sqlite3 [example.sqlite] |
| 終了  |  .exit <br> .quit|



## 専用コマンド

ドットから始まる専用のコマンドがあります。
[Command Line Shell For SQLite](https://sqlite.org/cli.html)


|   |   |
|---|---|
| 使用するDBファイルの指定<br>（ない場合は作成されます） |  .open example.sqlite |
| データベース一覧  |  .databases|
| テーブル一覧  |  .tables|
| テーブル情報  |  .schema table_name|


## 表のインポート

.modeを指定してからインポートします。タブ区切りなら`.mode tabs`、カンマ区切りなら`.mode csv`。指定の名前のテーブルがない場合は新規作成されます。
```
.mode tabs
.import example.tsv example_table
```

- テーブルがない場合（新規作成）：1行目がカラム名になり、データ型は全てtext型になります。
- テーブルが既にある場合：1行目もデータとして登録されます。カラム数が合わない場合エラーになります。データ型は違っても入るみたい（integerに文字列を入れても入った）。
- 空白や文字列"NULL"はtextになり、欠損値を入れることはできません。欠損値を入れたければインポート後に変更するとよいかも（e.g. UPDATE example_table SET col1 = NULL WHERE col1 = "";）


## データ型

[SQLiteで利用可能なデータ型](https://www.dbonline.jp/sqlite/type/index1.html)


## インデックス

[CREATE INDEX](https://www.sqlite.org/lang_createindex.html)


## pythonでの使用方法

組み込みライブラリsqlite3を呼び出すだけで使用できます。

DBへの接続/解除

```python
import sqlite3

# 接続（存在しない場合は新規作成されます）
conn = sqlite3.connect('example.sqlite')

# 接続解除
conn.close()
```

cursorオブジェクトを作って使う
```python
cur = conn.cursor()
cur.execute('SELECT * from example_table') # SQL送信
df = cur.fetchall() # 全ての行を取得
line = cur.fetchone() # 1行だけ取得

conn.close()
```

pandasで使う
```python
import sys
import pandas as pd

conn = sqlite3.connect('example.sqlite')

# 読み込み
df = pd.read_sql('SELECT * FROM example_table', conn)

# 書き込み
df.to_sql('table_name', conn, if_exists='replace', index=None)

conn.close()
```


















```bash
sqlite3 [example.sqlite]
```

終了

```sql
.exit or .quit
```



使用するDBファイルの指定（ない場合は作成されます）
```sql
.open example.sqlite
```

データベース一覧
```sql
.databases
```

テーブル一覧
```sql
.tables
```

テーブル情報
```sql
.schema table_name
```

テーブルのインデックス情報
```sql
.indices table_name
```

表示形式の設定
```sql
.mode list
.separator ", "
```


