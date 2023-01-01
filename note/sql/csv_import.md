# csvファイルのインポート(LOAD DATA INFILE)

[LOAD DATA INFILE](https://dev.mysql.com/doc/refman/5.6/ja/load-data.html)

ただしSQLのLOAD DATA INFILEはNULLになる条件などがややこしいので、pythonでpandasのto_sqlを使った方が分かりやすくておすすめ。




```sql
LOAD DATA [LOCAL] INFILE 'file.txt' -- mysqlサーバー上じゃなくてローカルにあるファイルの場合、`LOCAL`をつける。
  INTO TABLE table_name

  [FIELDS
    [TERMINATED BY ',']
    [OPTIONALLY ENCLOSED BY '"']
    ]

    [IGNORE 1 LINES]

    [(col_name_or_@user_var,...)] -- csvの列を左から順番に定義（DBのフィールドor変数に入力）
    [SET col_name = expr,...];
```
LOAD DATA INFILEの挙動
- 文字列型（Default: NULL）
  - 空白 -> 空白
  - NULL -> NULL（ちゃんとIS NULLで検出できるNULLになる）
  - NA   -> "NA"（文字列の"NA"。sqlにはNAがない）
- 数値型（Default: NULL）
  - 空白, NULL, NA -> 0（全部ゼロになる）

### SETの使い方

例）ファイルの最初のカラムを直接table_nameのcolumn1の値に使用し、2番目のカラムを除算演算の対象になるユーザー変数（@var1）に割り当てる。table_nameのcolumn2に@var1の演算結果を使用する。

```sql
LOAD DATA LOCAL INFILE 'file.txt'
  INTO TABLE table_name
  (column1, @var1)
  SET column2 = @var1/100;
```
例）入力値をダミーのユーザー変数に割り当て、その変数をテーブルカラムに代入しないようにしてその値を破棄する。

````sql
LOAD DATA LOCAL INFILE 'file.txt'
  INTO TABLE table_name
  (column1, @dummy, column2, @dummy, column3);
``````
例）NULLにする値を指定する

```sql
LOAD DATA LOCAL INFILE 'file.txt'
  INTO TABLE table_name
  (@var1, @var2, @var3)

  SET column1 = NULLIF(@var1, ''), # 空欄をNULL
      column2 = IF(@var2 = '' or @column = 'NA', NULL, @var2), # 空欄or'NA'をNULL
      column3 = NULLIF(@var3, 0) # 数値型の場合空欄も文字列も0に変換されるので、これでNULLにできる(0がないデータの場合)
  ```
