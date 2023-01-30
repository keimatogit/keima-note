# Linuxコマンドの結果をデータフレームとして読み込む



read.table（コマンド最後にパイプ(|)が必要）

```
df <- read.table("ls -l |", header = FALSE, sep = " ") %>%
	as_tibble()
```


read.csv（pipe()でコマンドを与える）

```
df <- read.csv(pipe("ls -l"), header = FALSE, sep = " ") %>%
	as_tibble()
```