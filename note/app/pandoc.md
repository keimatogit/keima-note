# pandoc

(221120) pandoc 2.19.2

```
-c conda-forge pandoc

# 確認
pandoc --version
```

markdownをhtmlに変換

[Pandocを使ってMarkdownを整形されたHTMLに変換する](https://qiita.com/cawpea/items/cea1243e106ababd15e7)

```
pandoc -s --toc -f markdown -t html -c styledheet.css -A footer.html index.md  >index.html
```