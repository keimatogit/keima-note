# Jupyter

```
conda create -n jupyter python=3.8 -y
conda activate jupyter
conda install -c conda-forge jupyter jupyterlab jupyter_contrib_nbextensions jupyterthemes pandas -y

# R関連
conda install -c conda-forge r=4.2 -y # 好きなバージョンを指定してください
conda install -c conda-forge r-base -y
conda install -c conda-forge r-essentials -y
conda install -c conda-forge r-irkernel r-openxlsx r-vegan -y
conda install -c bioconda bioconductor-biocinstaller bioconductor-qvalue r-argparse -y
```

Rを登録

```
R
IRkernel::installspec() # 今起動しているRを登録
IRkernel::installspec(name = 'ir42', displayname = 'R 4.2') # 好きな名前をつけて登録する場合
```


メモリ使用状況の表示(jupyter-resource-usage)
[Jupyter Notebookに現在のメモリ使用量を表示する](https://tarovlog.com/2021/04/10/jupyter-notebook-show-memory/)

```
conda install -c conda-forge jupyter-resource-usage
```


htimlにエクスポート
```
jupyter nbconvert --to html my-notebook.ipynb
# toc付き
jupyter nbconvert --to html_toc my-notebook.ipynb
```

```
conda install -c anaconda nbconvert=5.6.1=py36_0
```