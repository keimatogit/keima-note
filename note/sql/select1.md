# テーブルからデータを取得（SELECT）


[MySQL 5.6 リファレンスマニュアル 13.2.9 SELECT 構文](https://dev.mysql.com/doc/refman/5.6/ja/select.html)


<!-- @import "[TOC]" {cmd="toc" depthFrom=2 depthTo=6 orderedList=true} -->
<!-- code_chunk_output -->

1. [基本の形（SELECT ... FROM table_name）](#基本の形select-from-table_name)
2. [データ数（行数）を指定して取得（LIMIT）](#データ数行数を指定して取得limit)
3. [重複を除いて取得（DISTINCT）](#重複を除いて取得distinct)
4. [演算結果を取得（FUNCTION(colname)）](#演算結果を取得functioncolname)
5. [グループごとの演算結果の取得（GROUP BY）](#グループごとの演算結果の取得group-by)
6. [指定カラムや演算結果を任意の名前で取得（AS）](#指定カラムや演算結果を任意の名前で取得as)
7. [条件に合うデータの取得（WHERE）](#条件に合うデータの取得where)
8. [演算結果が条件に合うデータの取得（HAVING）](#演算結果が条件に合うデータの取得having)
9. [レコード（行）のソート（ORDER BY）](#レコード行のソートorder-by)
10. [DATE型のカラムから年や日付を抽出](#date型のカラムから年や日付を抽出)
11. [組み合わせた例](#組み合わせた例)
12. [SELECT結果をファイルに保存](#select結果をファイルに保存)

<!-- /code_chunk_output -->



[Sakila Sample Database](https://dev.mysql.com/doc/sakila/en/)を使用します。

employeesテーブル
```
+--------+------------+------------+-----------+--------+------------+
| emp_no | birth_date | first_name | last_name | gender | hire_date  |
+--------+------------+------------+-----------+--------+------------+
|  10001 | 1953-09-02 | Georgi     | Facello   | M      | 1986-06-26 |
|  10002 | 1964-06-02 | Bezalel    | Simmel    | F      | 1985-11-21 |
|  10003 | 1959-12-03 | Parto      | Bamford   | M      | 1986-08-28 |
+--------+------------+------------+-----------+--------+------------+
```

salariesテーブル
```
+--------+--------+------------+------------+
| emp_no | salary | from_date  | to_date    |
+--------+--------+------------+------------+
|  10001 |  60117 | 1986-06-26 | 1987-06-26 |
|  10001 |  62102 | 1987-06-26 | 1988-06-25 |
|  10001 |  66074 | 1988-06-25 | 1989-06-25 |
+--------+--------+------------+------------+
```

titlesテーブル
```
+--------+-----------------+------------+------------+
| emp_no | title           | from_date  | to_date    |
+--------+-----------------+------------+------------+
|  10001 | Senior Engineer | 1986-06-26 | 9999-01-01 |
|  10002 | Staff           | 1996-08-03 | 9999-01-01 |
|  10003 | Senior Engineer | 1995-12-03 | 9999-01-01 |
+--------+-----------------+------------+------------+
```


### 基本の形（SELECT ... FROM table_name）

指定カラムを取得
```sql
SELECT colname1, colname2, ... FROM table_name;
```

全てのカラムを取得
```sql
SELECT * FROM table_name;
```

文末のセミコロンの代わりに`\G`オプションをつけると、レコードごとの垂直表示になります
```sql
SELECT * FROM table_name \G
```

### データ数（行数）を指定して取得（LIMIT）

```sql
SELECT * FROM table_name LIMIT 10;  -- 10行だけ取得
```

### 重複を除いて取得（DISTINCT）

```sql
-- 指定カラムの重複を除く
SELECT DISTINCT colname FROM table_name;

-- 複数カラムの組み合わせの重複を除く
SELECT DISTINCT colname1, colname2, ... FROM table_name;

-- 全てのカラムの組み合わせの重複を除く
SELECT DISTINCT * FROM table_name;
```

### 演算結果を取得（FUNCTION(colname)）

|   |   |
|---|---|
| COUNT() | データ数（行数） |
| SUM() | 合計 |
| MAX() | 最大値 |
| MIN() | 最小値 |
| AVG() | 平均値 |
| STD() | 標準偏差 |
| VARIANCE() | 分散 |

COUNTは、カラム指定なし（COUNT(*)）の場合はNULLを含めた行数、カラム指定あり（COUNT(colname)）の場合はNULLを除いたデータ数を返します。

重複を除いたデータ数
```
```



### グループごとの演算結果の取得（GROUP BY）

```sql
SELECT FUNCTION(colname1) FROM table_name
 GROUP BY colname2, colname3, ...
```

```sql
SELECT COUNT(*) FROM table_name
　GROUP BY column_name1, column_name2;

SELECT AVG(column_name1) FROM table_name
　GROUP BY column_name2, column_name3;
```

### 指定カラムや演算結果を任意の名前で取得（AS）

```sql
SELECT colname1 AS my_colname1, AVG(colname2) AS col2_avg
 FROM table_name
 GROUP BY colname2;

SELECT AVG(colname1) AS col1_average FROM table_name;
```

### 条件に合うデータの取得（WHERE）
[SQLのwhereはinで複数条件を簡潔に記述できる サブクエリの併用も可](https://style.potepan.com/articles/23558.html)
```sql
SELECT * FROM table_name WHERE column_name = "character";
SELECT * FROM table_name WHERE column_name = num;
SELECT * FROM table_name WHERE column_name1 <= num1 OR column_name2 > num2;
SELECT * FROM table_name WHERE column_name BETWEEN start_value AND end_value; -- 境界値を含む
SELECT * FROM table_name WHERE column_name IS NULL;
SELECT * FROM table_name WHERE column_name IS NOT NULL;

-- ANDやORはINですっきり書ける
SELECT * FROM table_name WHERE column_name IN (num1, num2, num3);
SELECT * FROM table_name WHERE (column_name1, column_name2) IN ((num11, num21), (num12, num22), (num13, num23));
```

### 演算結果が条件に合うデータの取得（HAVING）

```sql
SELECT COUNT(*) FROM table_name GROUP BY column_name1 HAVING COUT(*) >= 100;
```

### レコード（行）のソート（ORDER BY）

`DESC`で降順。
```sql
SELECT * FROM table_name ORDER BY column_name;      
SELECT * FROM table_name ORDER BY column_name1, column_name2;
SELECT * FROM table_name ORDER BY column_name DESC;
SELECT * FROM table_name ORDER BY column_name1 DESC, column_name2;
```

### DATE型のカラムから年や日付を抽出
```sql
SELECT YEAR(date_column) FROM table_name; -- 年を取得
SELECT DATE_FORMAT(date_column, '%D') FROM table_name; -- 日付を取得
select DATE_FORMAT(date_column, '%Y-%m') FROM table_name; -- 指定フォーマットで取得
```


### 組み合わせた例
```sql
SELECT column_name1, column_name2, COUNT(*) AS N, AVG(column_name3 + column_name4) AS mean
FROM table_name
WHERE (column_name2 BETWEEN '2016-01-01' AND '2018-12-01') AND (column_name5 < 40)
GROUP BY column_name1, column_name2
ORDER BY mean DESC;
```

### SELECT結果をファイルに保存

[13.2.9.1 SELECT ... INTO 構文](https://dev.mysql.com/doc/refman/5.6/ja/select-into.html)
リダイレクトの方が楽かも。

```sql
SELECT a,b,a+b INTO OUTFILE '/tmp/result.txt'
  [FIELDS TERMINATED BY ','] [OPTIONALLY ENCLOSED BY '"']
  [LINES TERMINATED BY '\n']
  FROM table_name;
```
