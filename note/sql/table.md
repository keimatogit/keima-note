# テーブルの操作


<!-- @import "[TOC]" {cmd="toc" depthFrom=2 depthTo=6 orderedList=false} -->
<!-- code_chunk_output -->

- [テーブル情報の取得](#テーブル情報の取得)
- [テーブル作成（CREATE TABLE）](#テーブル作成create-table)
  - [設定を定義して作成](#設定を定義して作成)
  - [既存テーブルの設定をコピーして作成](#既存テーブルの設定をコピーして作成)
  - [既存テーブルの設定とデータをコピーして作成](#既存テーブルの設定とデータをコピーして作成)
- [テーブル情報の編集（ALTER TABLE）](#テーブル情報の編集alter-table)
- [テーブルのリセット（全レコードの削除）(TRUNCATE)](#テーブルのリセット全レコードの削除truncate)
- [テーブルの削除（DROP TABLE）](#テーブルの削除drop-table)
- [日付のデータ型のゼロについて](#日付のデータ型のゼロについて)

<!-- /code_chunk_output -->


## テーブル情報の取得

テーブル一覧の表示
```sql
show tables;
```

テーブルのフィールド情報取得
```sql
DESC table_name;
-- または
SHOW COLUMNS FROM table_name;
```
CREATE TABLE構文の形式でテーブル情報を表示
```sql
SHOW CREATE TABLE table_name;
```
レコード数を取得
```sql
SELECT COUNT(*) from table_name;
```

## テーブル作成（CREATE TABLE）

参考
[CREATE TABLE 構文](https://dev.mysql.com/doc/refman/5.6/ja/create-table.html)
[データ型](https://dev.mysql.com/doc/refman/5.6/ja/data-types.html)


### 設定を定義して作成

```sql
CREATE TABLE tabel_name (column_name data_type, ...)
```

例
```sql
DROP TABLE IF EXISTS user; -- すでに同じ名前のテーブルが存在している場合消しておく

CREATE TABLE user(

-- カラム名 データ型 NULL制約等
  id       INT(4) UNSIGNED ZEROFILL AUTO_INCREMENT,
  name     VARCHAR(10) NOT NULL,
  category TINYINT,
  birthday DATE, -- YYYY-MM-DD
  height   DECIMAL(3),
  weight   DECIMAL(4,1),

-- 主キーの指定（複数カラムの場合カンマ区切り。自動的にUNIQUE、NOT NULL制約がつき、インデックスが作成される）
  PRIMARY KEY (id),

-- 外部キー制約（自動的にインデックスが作成される）
  [CONSTRAINT [constraint_name]]
    FOREIGN KEY [index_name] (category) -- 外部キー制約をつけるカラム名
     REFERENCES parent_table(parent_table_column) -- 親テーブル名（親テーブルの対応カラム名）
     [ON UPDATE CASCADE] -- {RESTRICT,CASCADE,SET NULL}
     [ON DELETE CASCADE], -- {RESTRICT,CASCADE,SET NULL}

-- check制約
  [CONSTRAINT [constraint_name]]
    CHECK (height > 100)

) ENGINE=InnoDB DEFAULT CHARSET=utf8;
```

### 既存テーブルの設定をコピーして作成

```sql
CREATE TABLE newtable LIKE oldtable;
```

### 既存テーブルの設定とデータをコピーして作成
SELEC結果を使う。
```sql
CREATE TABLE newtable SELECT * FROM oldtable;
```


## テーブル情報の編集（ALTER TABLE）

テーブル名の変更
```sql
ALTER TABLE old_name RENAME TO new_name;
```

フィールド情報の変更
```sql
ALTER TABLE table_name MODIFY
  column_name data_type NOT NULL, ...
```

フィールドの追加と削除
```sql
ALTER TABLE table_name ADD COLUMN column_name data_type, ...;
ALTER TABLE table_name DROP COLUMN column_name;
```

主キーの指定と削除（カラムは残る）
```sql
ALTER TABLE table_name ADD PRIMARY KEY (column_name);
ALTER TABLE table_name DROP PRIMARY KEY;
```


## テーブルのリセット（全レコードの削除）(TRUNCATE)

```sql
TRUNCATE tabel_name; -- こちらの方が速い
-- または
DELETE FROM tabel_name;
```

## テーブルの削除（DROP TABLE）
```sql
DROP TABLE [IF EXISTS] tabel_name;
```
指定した名前のテーブルが存在しない場合、IF EXISTSなければエラーが出て、あればエラーが出ない。


## 日付のデータ型のゼロについて
[MySQL 5.6 リファレンスマニュアル  /  データ型  /日付と時間型](https://dev.mysql.com/doc/refman/5.6/ja/date-and-time-types.html)

>MySQLでは、DATE または DATETIME カラムに、日がゼロ、または月および日がゼロである日付の格納を許可しています。これは、正確な日付がわかっていない可能性のある生年月日を格納する必要があるアプリケーションで役立ちます。この場合は、単に日付を '2009-00-00' または '2009-01-00' として格納します。このような日付を格納する場合は、DATE_SUB() や DATE_ADD() などの完全な日付を必要とする関数で正しい結果が返されることは期待しないでください。

年・月のみのデータのときなど日付ゼロを入れられるが、pythonでデータをとってきた時に日付型だとエラーになってしまって扱いにくかった。ゼロをいれるのではなく、01などに統一しても良いと思う。
