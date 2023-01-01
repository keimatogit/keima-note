# カラム指定のgrep

awkを使う。
`awk -F "\t" 'パターン{アクション}'`


```
# 4列目が;Candida_albicansに一致する行を出力
cat 01_primary.tax | awk -F "\t" '$4==";Candida_albicans" {print $0}'

# 4列目に;Candida_albicansが含まれる行を出力
cat 01_primary.tax | awk -F "\t" '$4 ~ /;Candida_albicans/ {print $0}'

# 4列目に;Candida_albicansが含まれない行を出力
cat 01_primary.tax | awk -F "\t" '$4 !~ /;Candida_albicans/ {print $0}'

# 4列目と7列目に;Candida_albicansが含まれる行を出力
cat 01_primary.tax | awk -F "\t" '($4 ~ /;Candida_albicans/) && ($7 ~ /;Candida_albicans/) {print $0}'

# 4列目か7列目に;Candida_albicansが含まれる行を出力
cat 01_primary.tax | awk -F "\t" '($4 ~ /;Candida_albicans/) || ($7 ~ /;Candida_albicans/) {print $0}'
```

変数を使う
クオテーションがややこしいので注意
```
taxname="Candida_albicans"
cat 01_primary.tax | awk -F "\t" -v taxname=${taxname} '($4 ~ /;'"$taxname"'/) && ($7 ~ /;Candida_albicans/) {print $0}'
```