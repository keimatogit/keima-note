# publicデータのダウンロード（NCBI SRA/EMBL-EBI ENA）

## NCBI SRA (Sequence Read Archive database)

ラン番号（SRRxxxxxx）を指定してダウンロードする。

[sra-toolkit](https://github.com/ncbi/sra-tools/wiki)をインストールしておこう

```
conda create -n sra-tools install -c bioconda sra-tools
```

fastq-dump動作確認（ヘルプ表示）

```
fastq-dump --help
```

### プロジェクト番号からラン番号などのテーブルを取得

1. NCBIでBioProject （ex. PRJNA188481）を検索
2. SRA Experimentsをクリック
3. 全件表示して、Send to -> Run selecter
4. Select -> Download -> Metadataからメタデータのテーブルをダウンロード
5. カンマ区切り（しかもカラム名にカンマが含まれてる）たっだのでタブ区切りに変換`sed -e s/sea,region/sea_region/g SraRunTable.txt | sed -e s/,/\\t/g  > SraRunTable.tsv`

### ラン番号からfastqファイルを取得

標準出力で少しだけ見る（オプションのチェックなどに）

```
fastq-dump -X 5 -Z SRR769404
fastq-dump -X 5 -Z --defline-seq '@$ac $sn $sg $si $ri $rl' SRR769404
```

fastqのダウンロード

`--split-3`: ペアエンドを別ファイルに分ける
`--defline-seq`: 配列名のフォーマットを指定。'@$ac_$sn/$ri'だとR1は「SRR番号_連番/1」（R2は最後が/2）にする。（ex. デフォルト「@SRR769404.1 1 length=101」->「@SRR769404_5/2」）（221117: "fastq-dump" version 3.0.0 にて`sbatch ---wrap=""`の中だと、$acなどの正規表現が認識できないようで、「_/」だけとかになってしまった。シングルクオテーションとかクオテーションなしでもだめだった。wrap内でなければ大丈夫。）
`--defline-qual`: 3行目のフォーマットを指定。デフォは+シーケンス名。

```
fastq-dump \
  --gzip \
  --split-3 \
  --defline-seq '@$ac_$si/$ri' \
  --defline-qual '+' \
  --outdir fastq \
  SRR769404
```


## EMBL-EBI ENA

- [ENA](https://www.ebi.ac.uk/ena/browser/home)にアクセス
- 右上のEnter accessionにPRJEB38984を入れてViewをクリック
- サンプルメタデータ：Show Column Selection -> sample_title, library_strategyなど欲しい情報にチェックを入れる -> Download report TSV -> filereport_*.txtがダウンロードされる
- fastqファイル：ブラウザ上でチェックを入れてDownload selected files or Download Allでもよいし、↑のfastq_ftpやsubmitted_ftpのURLからwgetしてもよい

fastqはひとまずW0001-W0011の10人分をダウンロード(W0010は3時点しかなかったので除く)。

```
# メタデータからwget用のリストを作成
cat filereport_PRJEB38984.txt \
| sed '1d' \
| awk -F "\t" '{ if($14 ~ /^W000[123456789]_/ || $14 ~ /^W0011/ ) {print $12} }' \
| awk -F ";" '{print $1"\n"$2}' >wget_list.txt

# wgetでfastqをダウンロード（10人分だけでもめっちゃ時間かかる）
for i in $(cat wget_list.txt); do srun -p c6420 wget ${i}; done

# wget（WGSは特に時間がかかりすぎるのでsbatchで同時にwgetすることにした）
# でもあまりいっぺんにはできないみたい。。20ファイルずつくらいにする？
for i in $(cat wget_list.txt | grep "_WG_" | sed -n  "61, 80p")
do sbatch -c 8 -p c6420 -o wget_slurmout/%j.out --wrap="wget ${i}"
done
```
