# ファイルを結合するjoin

[JOINコマンドを使う際の注意](https://qiita.com/taruhachi/items/04ccca5e3bc20ad01791)

joinするにはキー列でsortしておく必要がある。sortとjoinで同じLANGを指定しておくと安心。

-tが区切り文字の指定。タブの指定は$'\t'になるようなので注意。

```
cat myfile.tsv | LANG=C sort -t$'\t' -k 1 >myfile_sort.tsv
```

```
LANG=C join -t$'\t' -j 1 \
  myfile_sort1.tsv \
  myfile_sort2.tsv
```


デフォルトではinner join（一致した行のみ出力）。`-a`オプションで一致しなかった行も出力させる。
- left join: `-a 1`
- right join: `-a 2`
- full (outer) join: `-a 1 -a 2`