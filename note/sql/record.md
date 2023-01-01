# レコード(行)の追加・編集・削除


<!-- @import "[TOC]" {cmd="toc" depthFrom=2 depthTo=6 orderedList=false} -->
<!-- code_chunk_output -->

- [追加（INSERT INTO）](#追加insert-into)
- [編集（UPDATE ... SET）](#編集update-set)
- [削除（DELETE FROM）](#削除delete-from)

<!-- /code_chunk_output -->


## 追加（INSERT INTO）

値を書いて追加
全列に入れる場合はカラム名を省略できる。文字列はクオテーションで囲む。
```sql
INSERT INTO table_name [(column1, column2, ...)] VALUES (value1, value2, ...);
```

別テーブルのSELECT結果を追加
```sql
INSERT INTO table_name [(column1, column2, ...)]
  SELECT column1_2, column2_2, ...
  FROM table_name2
  WHERE column1_2 > 100;
```

## 編集（UPDATE ... SET）

全部の行の値を変更
```sql
UPDATE table_name SET column1=1000, column2=10;
```

条件に合う行の値を変更
```sql
UPDATE table_name SET column1=1000, column2=10 WHERE id=10;
```

他の列の演算結果を格納
```sql
UPDATE table_name SET column2=(column1*0.1);
```

date型のカラムに日付を１日足す
```sql
UPDATE table_name SET birthday = (birthday + INTERVAL 1 DAY);
```

## 削除（DELETE FROM）

条件に合う行を削除
```sql
DELETE FROM table_name WHERE column1 >= 100 AND column2 >= 10;
```

全ての行を削除
```sql
DELETE FROM table_name;
-- または
TRUNCATE TABLE; -- こっちの方が速い
```
