# MySQLサンプルデータの使い方

こちらからダウンロードできます：[MySQLサンプルデータ](https://dev.mysql.com/doc/index-other.html)

解凍して読み込みます
```shell
wget https://downloads.mysql.com/docs/world.sql.gz
gunzip world.sql.gz
mysql -u root -p < world.sql
```