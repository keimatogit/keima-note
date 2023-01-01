# データベースの操作

MariaDBやMySQLのDBサーバーには、名前がついた「データベース」が複数あり、各「データベース」に複数のテーブルが格納されている構造です。DBを使用する際はどの「データベース」を使用するか、データベース名を指定します。

DB一覧の表示
```sql
SHOW DATABASES;
```

使用DBの指定
```sql
USE database_name;
```

DBの作成
（すでにその名前のDBがある場合は何も起きない。IF NOT EXISTSがなければエラーが出て、あればエラーが出ない。）
```sql
CREATE DATABASE [IF NOT EXISTS] database_name;
```

DBの削除
（その名前のDBがない場合、IF EXISTSがなければエラーが出て、あればエラーが出ない。）
```sql
DROP DATABASE [IF EXISTS] database_name;
```
