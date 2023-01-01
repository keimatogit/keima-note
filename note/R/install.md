

R-baseとRToolsをインストール(https://cran.ism.ac.jp/)？

linux
> Download and install the latest base R version via the Linux package manager (r-base which consists of r-base-core and r-base-devel) 

## ubuntuに最新のRをインストール

[Ubuntu20.04にRをインストールする方法](https://www.trifields.jp/how-to-install-r-on-ubuntu-2004-4335)

```
sudo apt update # パッケージの一覧を更新
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 # GPGキーを登録
sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" # aptのリポジトリに登録
sudo apt install --no-install-recommends r-base
```