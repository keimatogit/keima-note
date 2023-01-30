# EggNOG-mapper

> EggNOG-mapper は、新規配列の機能アノテーションを高速に行うツールである。・・・eggNOG-mapperの一般的な用途としては、新規ゲノム、トランスクリプトーム、あるいはメタゲノム遺伝子カタログのアノテーションが挙げられます。

基本的にDIAMOND blastpで探すみたい（blastxもできる）。

[github](https://github.com/eggnogdb/eggnog-mapper)
[github - eggNOG mapper v2.1.5 to v2.1.8](https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.8)


## セットアップ

eggnog-mapperのインストール

```
conda create -n eggnog-mapper -c bioconda eggnog-mapper

# 確認
conda activat eggnog-mapper
emapper.py -h 
```

`download_eggnog_data.py`でレファレンスデータベースをダウンロードします。eggNOG、eggNOGのtaxa、Diamond blast用のeggNOGタンパク質データベースがダウンロードされる。毎回[y/n]を聞かれるので、eggNOGダウンロードにすごく時間かかった後にtaxaとdiamondについてもきかれるんん注意。（230125: srunでかけると止まってしまって進まなかった。初めにy/nを聞かれるからその部分でうまく進まなかったのかも。srunなしでそのまま実行したら大丈夫だった。）

```
# オプションを確認しよう
download_eggnog_data.py -h

# ~/db/eggnog_mapperにダウンロード 
download_eggnog_data.py \
 --data_dir ~/db/eggnog_mapper
 
# Pfamとかもダウンロードする場合はオプションをつけるよう言われます
## Skipping novel families diamond database (or already present). Use -F and -f to force download
## Skipping Pfam database (or already present). Use -P and -f to force download
## Skipping MMseqs2 database (or already present). Use -M and -f to force download
## No HMMER database requested. Use "-H -d taxid" to download the hmmer database for taxid
## Finished.
```

４つのファイルがダウンロードされました。
- eggnog.db
- eggnog_proteins.dmnd
- eggnog.taxa.db
- eggnog.taxa.db.traverse.pkl


（補足）データベースの場所について： 環境変数`EGGNOG_DATA_DIR`をセットしておけば、eggnog-mapperの全てのコマンドで`--data_dir`を指定する必要がなくなります（デフォルトがEGGNOG_DATA_DIRになる）。たとえば`/home/user/db/eggnog_mapper`を`EGGNOG_DATA_DIR`にセットする場合は、~/.bashrcに以下のように書き加えます。EGGNOG_DATA_DIRがセットされていない場合は、eggnog-mapperディレクトリが含まれているフォルダ内の`data`フォルダがデフォルト値になります（condaでggnog-mapperをインストールした場合、`miniconda3/envs/eggnog-mapper/lib/python3.10/site-packages/data/eggnog.db`などを探しに行く）。

```
export EGGNOG_DATA_DIR=/home/user/db/eggnog_mapper
```