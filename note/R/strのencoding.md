# strのencoding

システムのlocaleによる。変更してよければlocaleを変更してしまうといいのかも。

localeの確認
```
Sys.getlocale()
```

winの932だとEncodingがunkwonになるみたい。macやubuntuだとutf-8となる。
Encodingがunknownの場合にShit-jisからutf-8に変えてurlencodeにする例
```
companye <- "ユニクロ"
if (Encoding(company_name) == "unknown"){
  company_name <- iconv(company_name,'SHIFT_JIS','UTF-8')
}
company_name_urlencode <- URLencode(company_name)
```

withrパッケージでsortの順番などは変えられる。
withrでEncodingをutf-8にする方法はあるのかな？色々やってみたけどだめだった。
```
tmp <- c("あ", "イ", "う", "エ", "お")
sort(tmp)
withr::with_collate("C", sort(company_name))
```