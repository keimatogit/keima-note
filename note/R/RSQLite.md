# RSQLite

Rでsqliteを使う

```r
library(RSQLite)

con <- dbConnect(SQLite(), db, synchronous="off")
sql <- paste0("SELECT gene_id, knumber from dat WHERE gene_id IN ('", str_flatten(unique(genes), "','"), "')")
sql_result <- dbGetQuery(con, sql)
```