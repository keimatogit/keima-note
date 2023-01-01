# awk内でbash変数を使う

正直よくわからないので注意。
都度やってみてください。

```
taxname="Calbicans"
echo | awk '{print "'$taxname'"}'
```

```
cat ${sample}_primary.tax \
 | awk -F "\t" '($4 ~ /;'"$tax2"'/) && ($7 ~ /;'"$tax2"'/) {print $1}' 
```