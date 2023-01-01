# jupyter-notebookでRを使う


## jupyter-notebookでRを使う
[jupyter notebookでRを使う](https://nxdataka.netlify.app/rjup/)

jupyter用のRカーネルパッケージをインストールして登録

Rで実行
```r
install.packages('IRkernel')
IRkernel::installspec() # 登録
IRkernel::installspec(name = 'ir42', displayname = 'R 4.2') # 好きな名前をつけて登録する場合
```

condaの場合
```
conda install -c r r r-irkernel
```

## プロットサイズ

jupyterのR環境はreprパッケージで管理されており、reprのオプション設定でプロット表示サイズを調整できる。

```r
options(repr.plot.width=7, repr.plot.height=7, repr.plot.res=120) # デフォルト値
```

## ggplotでの日本語表示

ggplotで日本語を使用すると警告がたくさん出るので、`suppressWarnings(plot(g))`のように非表示にすると良い。ちなみにggplotで日本語フォントを使うには`theme(text = element_text(family = "HiraginoSans-W3"`のようにtheme内で指定する。

使えるフォントは`systemfonts::system_fonts()`で確認
良いのがない場合はaptで追加
```
sudo apt install fonts-ricty-diminished
```
