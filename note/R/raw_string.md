# raw string

[R 4.0.0 にさっさと移行しようと思ったたった一つの理由：新しいraw string記法【 r"" 】](https://qiita.com/taiyodayo/items/250483de6228eb298c80)
[Pythonでエスケープシーケンスを無視（無効化）するraw文字列](https://note.nkmk.me/python-raw-string-escape/)

""の代わりにr"()"で文字列を囲むと、特殊文字をエスケープなしで記入できます（自動的にエスケープを入れてくれる）。

`awk -F "\t" '{print $1"\t"$2}'`を記入してみましょう。

```
example <- r"(awk -F "\t" '{print $2"\t"$6}')"

example
## [1] "awk -F \"\\t\" '{print $2\"\\t\"$6}'"

cat(example)
## awk -F "\t" '{print $2"\t"$6}'
```