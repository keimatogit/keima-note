# miniconda環境

## miniconda3のインストール

[minicondaのサイト](https://docs.conda.io/en/latest/miniconda.html)から最新を取得

```shell
# wget
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
# getがないとき
curl -O https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

bash Miniconda3-latest-Linux-x86_64.sh
# インストールが始まるのでyesとか言う
```

## チャンネル登録と確認
```shell
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda

# 確認
conda config --get channels
```

## 仮想環境の作成など

仮想環境を作成

```shell
conda create -n myenv [python=3.8]
```

仮想環境の有効化

```
conda activate myenv
```

仮想環境を終了

```
conda deactivate
```

仮想環境を削除

```
conda remove -n myenv --all
```

仮想環境の一覧

```
conda info --envs
```

## ソフトウェアの検索とインストール

`-y`で途中の質問を省略
```
conda search software_name
conda install software_name [-y]
conda install tensorflow=1.13.1 # バージョン指定
conda install tensorflow=1.13.1=build # バージョン＆build指定
```

## minicondaのアンインストール
```
rm -rf ~/miniconda
rm -rf ~/.condarc ~/.conda ~/.continuum
# さらに~/.bash_profileを編集
```


### Jupyter notebbokと拡張機能

新しいnbconvertではdownload html with tocができず、5.6.1にしたらできた。でもdownload html with tocだとfigが別ファイルになるので、普通のdownload htmlの方がいいかも
```
conda install jupyter jupyter_contrib_nbextensions -y
conda install nbconvert=5.6.1
conda install samtools=1.14=hb421002_0
```

### pythonライブラリ
```
# データベース関連
conda install mysql-connector-python sqlalchemy -y
# ほか
conda install pandas seaborn ipywidgets -y
```

### R関連

[jupyter notebookでRを使う](https://nxdataka.netlify.app/rjup/)

パッケージのインストールは`conda install`とR内での`install.packages("")`を併用していたらヘンになったので、どちらかに統一した方がよさそう。

```
# R本体と基本パッケージ
conda install -c conda-forge r r-essentials r-irkernel  # -c r のは古かったので

# 追加パッケージ例
conda install -c conda-forge r-argparser # 引数受け取り
conda install -c conda-forge r-vegan r-coin  r-effsize # 統計とか
conda install -c conda-forge r-googlesheets4 # google sheetを扱う
conda install -c bioconda bioconductor-limma # limma

# DB関連パッケージ（SQLite, MySQL, MariaDB）
#（rmysqlやrmariadbはMySQLやMariaDBが入ってないとインストールできないみたい）
# r-mariadbはr-dbiも必要だけどr-essentialで既に入ってるはず）
conda install -c conda-forge r-rsqlite r-rmysql r-rmariadb

# rstudioを入れる場合
conda install -c r rstudio
```

