---
title: "Ubuntu setup memo"
author: "matsumoto"
date: 2021/02/04
output: html
html:
  toc: true
---

# WindowsでUbuntu環境構築

<!-- @import "[TOC]" {cmd="toc" depthFrom=2 depthTo=6 orderedList=false} -->
<!-- code_chunk_output -->

- [WSL2とUbuntuのインストール](#wsl2とubuntuのインストール)
- [Windows terminalを使う](#windows-terminalを使う)
- [aptコマンドでのソフト管理](#aptコマンドでのソフト管理)
- [マウント](#マウント)
  - [サーバーのマウント](#サーバーのマウント)
  - [Googleドライブのマウント](#googleドライブのマウント)
- [環境構築](#環境構築)
  - [condaを使う（おすすめ）](#condaを使うおすすめ)
  - [condaを使わない場合](#condaを使わない場合)
    - [R関連](#r関連)
    - [MariaDB](#mariadb)
    - [python関連](#python関連)

<!-- /code_chunk_output -->


## WSL2とUbuntuのインストール
[Windows 10 用 Windows Subsystem for Linux のインストール ガイド](https://docs.microsoft.com/ja-jp/windows/wsl/wsl2-kernel)

1. プラグラムと機能 -> Windowsの機能の有効化または無効化 で「Linux用Windowsサブシステム」「仮想マシン プラットフォーム」を有効にしてPCを再起動。
2. [Windows Subsystem for Linux インストール ガイド 手順4](https://docs.microsoft.com/ja-jp/windows/wsl/wsl2-kernel)のリンクから「WSL2 Linux カーネル更新プログラム パッケージ」をダウンロード、起動してインストール。
3. Power Shellでコマンド `wsl --set-default-version 2`を実行し、WSL2をデフォルトのバージョンに指定。
4. Microsoft storeからUbuntu(LTSの最新)をインストール。起動してユーザー名とパスワードを設定します。

WSLのバージョン確認や設定はWindows Powershellでいつでもできます
```
wsl --list --verbose             # 現在のバージョン確認
wsl --set-default-version 2      # WSL2をデフォルトに設定
wsl --set-version Ubuntu-20.04 2 # 現在のubuntuのバージョンをWSL2に変更
```

## Windows terminalを使う
[Windowsで至高のターミナル生活を求めて(Windows Terminal編)](https://www.asobou.co.jp/blog/web/windows-terminal)

Windows terminalでubuntuを使うとタブ分割ができて便利です。見た目のカスタマイズもできます。

Microsoft storeからWindows terminalをインストールします。起動して画面上で Ctrl + , と入力すると設定画面が開きます。「JSONファイルを開く」で直接settings.jsonファイルを開いて好きな設定を書き込むこともできます。

settings.jsonの書き方
- デフォルトのプロファイルは`"defaultProfile": `で指定します。`"profiles":`内の`"list":`のなかからデフォルトにしたいプロファイルを探して、その`"guid": `を指定します。例： `"defaultProfile": "{07b52e3e-de2c-5db4-bd2d-ba144ed6c273}",`
- デフォルトの見た目の設定は`"profiles":`内の`"defaults":` の中に追加します。例：`"fontSize": 14, "colorScheme": "Tango Dark"`
- wslでubuntuを起動した時の場所をubuntuのホームディレクトリにするには、`"profiles":`内の`"list":`のなかのubuntuの欄に`"startingDirectory": `を追加して指定します。例：`"startingDirectory" : "//wsl$/Ubuntu-20.04/home/user_name"`



## aptコマンドでのソフト管理

[aptコマンドチートシート](https://qiita.com/SUZUKI_Masaya/items/1fd9489e631c78e5b007)

はじめにパッケージ一覧の更新と入っているパッケージのアップデートをしておきます
```
sudo apt update # パッケージ一覧の更新
sudo apt upgrade # パッケージのアップデート
```

ソフトをインストールするとき
```
sudo apt update # インストール前にパッケージ一覧を更新
sudo apt search package_name # パッケージを検索(部分一致)
sudo apt install package_name1 package_name2 package_name3 # インストール（複数可）
```

アップデートするとき
```
sudo apt upgrade
```

アンインストールするとき
```
sudo apt remove package_name
sudo apt --purge remove package_name # 依存関係があるパッケージを含めて完全に削除
```

## マウント

ローカルドライブは `/mnt/`  に自動的にマウントされます。

### サーバーのマウント

```
# マウントポイントを作成しておく
sudo mkdir /mnt/mount_folder

# マウント（-o ro：読み込み専用でマウント）
sudo mount [-o ro] -t drvfs '\\192.168.1.xx\xxx\xxx' /mnt/mount_folder

# マウント解除（必要なら）
sudo umount /mnt/mount_folder
```


### Googleドライブのマウント

google-drive-ocamlfuseを使います。

[google-drive-ocamlfuseでGoogleドライブをマウント](https://qiita.com/cabbage_lettuce/items/c4544b3e5cd28caf04bc)
```
# インストール
sudo add-apt-repository ppa:alessandro-strada/ppa
sudo apt update
sudo apt install google-drive-ocamlfuse

# 認証（アクセス用の情報が~/.gdfuse/default/に保存される）
google-drive-ocamlfuse

# マウントポイントを作ってマウント
mkdir ~/GoogleDrive
google-drive-ocamlfuse ~/GoogleDrive

# マウント解除（必要なら）
fusermount -u ~/GoogleDrive
```

GUIブラウザがないときの認証方法：
[Headless Usage & Authorization](https://github.com/astrada/google-drive-ocamlfuse/wiki/Headless-Usage-&-Authorization)
[sshごしにgoogle-drive-ocamlfuseのOAuth認証を行う方法](http://moguno.hatenablog.jp/entry/2016/03/24/010502)

Windows terminalではできなかったので、Ubuntu自体を起動して実行するとよいかも。
```
echo $'#!/bin/sh\necho $* > /dev/stderr' > xdg-open
chmod 755 xdg-open
PATH=`pwd`:$PATH google-drive-ocamlfuse
```
表示されたURLにアクセスして認証を行います。しばらく待ってubuntuでAccess token retrieved correctly.と表示されたら認証終了です。

## 環境構築

 python関連はconda、Rは色々なパッケージを使う場合はconda外で管理するのがおすすめ（conda標準セット（r-essentials）意外に色々なパッケージを入れていたらしばしば不具合が出たので）。


### condaを使う（おすすめ）

[condaの使い方](conda.md)

### condaを使わない場合

#### R関連

- r-base-core：R本体
- r-recommended：CRANのRパッケージに依存するメタパッケージ
- r-base：r-base + r-recommended：Rstudioにはこれが必要？

```
sudo apt update
sudo apt install r-base-core

# パッケージのインストールに以下の３つが必要
sudo apt install libcurl4-openssl-dev
sudo apt install libssl-dev
sudo apt install libxml2-dev

# RMariaDBのインストールに必要みたい
sudo apt-get install libmariadb-dev

# systemfontsのインストールに必要みたい
sudo apt install libfontconfig1-dev
```

Rを起動して必要なRパッケージのインストール
```r
install.packages("tidyverse")
install.packages("RMariaDB")
install.packages("readxl")
install.packages("googlesheets4")

# jupyter-notebookでRを使う場合
install.packages('IRkernel')
IRkernel::installspec()

# ggplotで日本語使う場合あると便利
install.packages('systemfonts')
```

Rstudio Server：
R開発環境。起動してからwindowsのブラウザで `http://localhost:8787` にアクセスすると使用できる。このとき聞かれるユーザー名とパスワードはubuntuのもの。
[Download RStudio Server for Debian & Ubuntu](https://rstudio.com/products/rstudio/download-server/debian-ubuntu/)

```
# 必要なソフトをインストール
sudo apt update
sudo apt install r-base
sudo apt install gdebi-core

# RStudio Serverをダウンロード（バージョン名は新しいものを入れること）
wget https://download2.rstudio.org/server/bionic/amd64/rstudio-server-1.4.1103-amd64.deb

# インストール
sudo gdebi rstudio-server-1.4.1103-amd64.deb

# 起動
sudo rstudio-server start
```


#### MariaDB
[MariaDB Package Repository Setup and Usage](https://mariadb.com/kb/en/mariadb-package-repository-setup-and-usage/)
```
sudo apt install apt-transport-https
curl -LsS https://downloads.mariadb.com/MariaDB/mariadb_repo_setup | sudo bash

sudo apt install mariadb-server
sudo mysql_secure_installation
```

#### python関連

python関連はconda使用がおすすめ！

pip：
Pythonパッケージのインストールなどを行う
```
sudo apt install python3-pip

# pip3 install使用時にパスを通すよう言われるかも。
echo export PATH=$PATH:/home/user_name/.local/bin >> ~/.bashrc
```

jupyter notebook：
Pythonなどの開発環境
```
sudo apt update
sudo apt install jupyter-notebook  
jupyter --version # 確認
jupyter-notebook --no-browser # 起動。表示されるURLをWindowsのウェブブラウザに貼り付けて使う。
```

ToC2：
jupyter notebook拡張機能。サイドバーにTOCを表示できる
```
pip3 install jupyter_contrib_nbextensions
jupyter contrib nbextension install --user
pip3 install jupyter_nbextensions_configurator
jupyter nbextensions_configurator enable --user
```

jupyterthemes：
jupyter notebookデザイン調整。フォントファミリーはchromの設定->デザイン->フォントをカスタマイズ->固定幅フォント欄で設定できる。
[7 Essential Tips for Writing With Jupyter Notebook](https://towardsdatascience.com/7-essential-tips-for-writing-with-jupyter-notebook-60972a1a8901)
```
pip3 install jupyterthemes  
pip3 install --upgrade jupyterthemes

# 好きなテーマやフォントサイズを指定（-N, -Tはノート名とツールバーの表示を指定）
jt -t oceans16 -ofs 11 -T -N
jt -fs 13 -ofs 13 -f "Ricty Diminished" -T -N

# テーマをリセット
jt -r
```
-tf text/markdown cell font
-fs code font size
-ofs output area font size


ライブラリ：
必要な物をpip3でインストールする。

```
pip3 install pandas numpy sqlalchemy mysql-connector-python
```

mariadbパッケージを入れるにはaptでlibmariadb-devを入れておく必要があるみたい
（参考：[Problem with pip install mariadb - mariadb_config not found](https://stackoverflow.com/questions/63027020/problem-with-pip-install-mariadb-mariadb-config-not-found)）
```
sudo apt update -y
sudo apt install libmariadb-dev
```
