# awk内でbash変数を使う

```
# 14列目が指定のパターンを含む場合12列目をプリント
cat filereport_PRJEB38984.txt | sed '1d' | awk -F "\t" '{ if($14 ~ /^W000[123456789]_/ || $14 ~ /^W0011/ ) {print $12} }' 
```