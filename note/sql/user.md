
# DBユーザーの管理

DBサーバーのユーザー管理です。各ユーザーには、データベース・テーブルごとに権限をつけられます。
ユーザー名とホスト名はそれぞれシングルクオテーションで囲みます（ 'user01'@'hostname'など）。ホスト名'%'は全てのホストにマッチします（どんなホストでも接続できる）。ユーザー名のみのアカウントは'user_name'@'%'と同等に扱われます。

<!-- @import "[TOC]" {cmd="toc" depthFrom=2 depthTo=6 orderedList=false} -->
<!-- code_chunk_output -->

- [ユーザー一覧と権限の確認](#ユーザー一覧と権限の確認)
- [ユーザーの追加・削除・変更](#ユーザーの追加削除変更)
- [ユーザー権限の変更](#ユーザー権限の変更)
  - [権限の付与](#権限の付与)
  - [権限の削除](#権限の削除)

<!-- /code_chunk_output -->


## ユーザー一覧と権限の確認

ユーザー一覧
```sql
SELECT host,user,password FROM mysql.user;
```

ログインユーザーの確認
```sql
SELECT USER();
```

ユーザー権限確認
```sql
SHOW GRANTS for 'user01'@'%';
```


## ユーザーの追加・削除・変更

ユーザーの追加と削除（パスワード設定は任意）
```sql
CREATE USER 'user01'@'%' [IDENTIFIED BY 'mypass'];
DROP USER 'user01'@'hostname';
```

ユーザー名の変更

```sql
RENAME USER old_name TO new_name;
```

パスワード変更
（"FOR ユーザー名"を省略した場合、ログインユーザーのパスワードを変更）
```sql
SET PASSWORD [FOR user01] = PASSWORD('mypass');
```


## ユーザー権限の変更

権限の付与・削除は、GRANT OPTIONがついてるユーザーから、そのユーザーと同等以下の権限のみ可能

### 権限の付与

[MySQL 5.6 リファレンスマニュアル  /  ...  /  GRANT 構文](https://dev.mysql.com/doc/refman/5.6/ja/grant.html)

```sql
GRANT [権限] ON [DB名].[テーブル名] TO '[ユーザー名]'@'[ホスト名]';
```
```sql
# 特定のテーブルでのINSERT権限を与える
GRANT INSERT ON database_name.table_name TO 'user01'@'localhost';

# 全てのテーブルで読み込み権限のみ与える
GRANT SELECT ON *.* TO 'user01'@'localhost';

# 全ての権限をつける(GRANT OPTION を除く)
GRANT ALL ON *.* TO 'user01'@'localhost';

# 最後にWITH GRANT OPTIONをつけるとほかのアカウントへの権限付与・削除もできるようになる
GRANT ALL ON *.* TO 'user01'@'localhost' WITH GRANT OPTION;
```

### 権限の削除

[MySQL 5.6 リファレンスマニュアル  /  ...  /  REVOKE 構文](https://dev.mysql.com/doc/refman/5.6/ja/revoke.html)

```sql
REVOKE [権限] ON [DB名].[テーブル名] FROM '[ユーザー名]'@'[ホスト名]';
```
```sql
# INSERT権限を削除
REVOKE INSERT ON *.* FROM 'testuser'@'%';

# GRANT OPTIONを削除
REVOKE GRANT OPTION ON *.* FROM testuser@'%';

# 全ての権限とGRANT OPTIONを削除
REVOKE ALL PRIVILEGES, GRANT OPTION FROM 'testuser'@'%';
```
