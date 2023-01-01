# インデックス(key)について

[MySQL のインデックスの使用の仕組み](https://dev.mysql.com/doc/refman/5.6/ja/mysql-indexes.html)

- インデックス = キー
- インデックスを作成したカラムは値の検索がすばやくできる。
- Primary Key(主キー)に指定したカラム、Foreign Key(外部キー)制約を指定したカラムには、自動的にインデックスが作成されるので、明示的にインデックスを作成する必要はない。


<!-- @import "[TOC]" {cmd="toc" depthFrom=2 depthTo=6 orderedList=false} -->
<!-- code_chunk_output -->

- [インデックス情報の取得](#インデックス情報の取得)
- [インデックスの追加（ALTER TEBLE ... ADD INDEX）](#インデックスの追加alter-teble-add-index)
- [インデックスの削除（ALTER TEBLE ... DROP INDEX）](#インデックスの削除alter-teble-drop-index)

<!-- /code_chunk_output -->



## インデックス情報の取得

テーブル中の詳しいインデックス情報を取得（Key_nameがインデックス名）
[SHOW INDEX 構文](https://dev.mysql.com/doc/refman/5.6/ja/show-index.html)


```sql
SHOW INDEX FROM table_name;
```

データベース中のインデックス情報を取得

```sql
SELECT TABLE_SCHEMA, TABLE_NAME, COLUMN_NAME, INDEX_NAME
FROM INFORMATION_SCHEMA.STATISTICS
WHERE TABLE_SCHEMA='database_name';

# テーブルを限定して
SELECT TABLE_SCHEMA, TABLE_NAME, COLUMN_NAME, INDEX_NAME
FROM INFORMATION_SCHEMA.STATISTICS
WHERE TABLE_NAME='table_name';
```

## インデックスの追加（ALTER TEBLE ... ADD INDEX）
複合インデックスの場合はカラム名をカンマ区切りで指定。

```sql
ALTER TABLE table_name ADD INDEX [index_name](column_name);
ALTER TABLE table_name ADD INDEX [index_name](column_name1, column_name2);
```

## インデックスの削除（ALTER TEBLE ... DROP INDEX）
```sql
ALTER TABLE table_name DROP INDEX index_name;
```
