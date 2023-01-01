# データベースを使用する



## RMariaDBパッケージを使う
[RMariaDB](https://rmariadb.r-dbi.org/)

```r:rmariadb.R
library(RMariaDB)
library(DBI)

# connect to DB
conn <- dbConnect(MariaDB(),
                  user       = "user_name",
                  password   = "password",
                  dbname     = "db_name",
                  host       = "host_name")

# テーブル全体を取得
df_all <- dbReadTable(conn, "table_name")

# フィールド名を取得
field_names <- dbListFields(conn, "table_name")

# SQL実行結果を取得
query <- "SELECT * FROM table_name LIMIT 5"
res <- dbSendQuery(con, query)
df <- dbFetch(res)
dbClearResult(res)

# データフレームをDBへインポート
dbWriteTable(con,
             value     = df,
             name      = "table_name",
             row.names = FALSE,
             append    = TRUE,
             overwrite = FALSE)

# 最後に接続解除
dbDisconnect(con)
```

## RMySQLパッケージを使う

```r
library(DBI)
library(RMySQL)
```
