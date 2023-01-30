# MetaPhlAn/StrainPhlAn

[GitHub](https://github.com/biobakery/MetaPhlAn)
[GitHub - MetaPhlAn 4](https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-4)
[GitHub - StrainPhlAn 4](https://github.com/biobakery/MetaPhlAn/wiki/StrainPhlAn-4)

```
conda create --name metaphlan4 -c conda-forge -c bioconda python=3.7 metaphlan
```

MetaPhlAn4のデータベースはconda外のフォルダに置くことが推奨されているらしい（大きいから？）。配列などをダウンロードした後bowtie2のインデックスを作製するので時間がかかる。

> データベースを別の場所にインストールした場合は、-bowtie2db <データベースフォルダ>を使用してMetaPhlAnを実行することを忘れないでください!

```
# 最新のMetaPhlAnDBをダウンロード
metaphlan --install --bowtie2db <database folder>
```

上のコマンドだとmarker_info.txtやspecies.txtがダウンロードされなかったので、別途[Segatalab FTP](http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/)がらダウンロードした。

```
wget http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vJan21_CHOCOPhlAnSGB_202103_marker_info.txt.bz2
wget http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vJan21_CHOCOPhlAnSGB_202103_species.txt.bz2
bunzip2 mpa_vJan21_CHOCOPhlAnSGB_202103_marker_info.txt.bz2
bunzip2 mpa_vJan21_CHOCOPhlAnSGB_202103_species.txt.bz2
```

ちなみに含まれているpickleファイル（mpa_vJan21_CHOCOPhlAnSGB_202103.pklpython）をpythonで読み込むには以下。中身は入れ子になったdict。

```
import os
import pickle
import bz2
from Bio import SeqIO
db = pickle.load(bz2.BZ2File("mpa_vJan21_CHOCOPhlAnSGB_202103.pkl"))
db.keys()
## dict_keys(['taxonomy', 'markers', 'merged_taxon'])

# db['markers']がmarker_info.txtに相当してるみたい
```

(2022-11-11)~/db/metaphlan/bowtie2db/221111にダウンロード。

conda内のデフォルトの置き場所は、miniconda3/envs/metaphlan4/lib/python3.7/site-packages/metaphlan/metaphlan_databases。metaphlan実行時に`-bowtie2db`でdb場所を指定しなかった場合、この場所に自動的にダウンロード・インデックス作製される。



## StrainPhlAn

事前にアダプター除去とクオリティコントロールを行なっておく。

MetaPhlAn: marker database of MetaPhlAnにマッピングしてsamファイルを作成

```
mkdir sams profiles bowtie2
metaphlan \
  ${sample}_R1.fastq.gz,${sample}_R2.fastq.gz \
  --input_type fastq \
  --nproc 8 \
  --bowtie2db ~/db/metaphlan4/bowtie2db \
  [--unclassified_estimation \]
  -s sams/${sample}.sam.bz2 \
  -o profiles/${sample}_profiled.tsv \
  --bowtie2out bowtie2/${sample}.bowtie2.bz2

# 以前にMetaPhlAnを実行してbowtie2outファイルがある場合は、すぐに結果が出る
metaphlan \
  metagenome.bowtie2.bz2 \
  --input_type bowtie2out \
  --nproc 8 \
  -o profiled.txt

# 独立にbowtie2を実行してsamファイルから読み込むこともできる
bowtie2 \
  --sam-no-hd \
  --sam-no-sq \
  --no-unal \
  --very-sensitive \
  -S metagenome.sam \
  -x metaphlan_databases/mpa_vJan21_CHOCOPhlAnSGB_202103 \
  -U metagenome.fastq
metaphlan metagenome.sam \
  --input_type sam \
  -o profiled_metagenome.txt
```

sample2markers.py: StrainPhlAnの入力となるコンセンサスマーカーファイルを作成。正規表現をつかって複数サンプルのsamをいっぺんにかけられる。
(2022-11-24)時間がかかっているので１サンプルずつで別々にsbatchした方がいいかも！

> 各サンプルについて、その中に含まれる全種類の菌株を再構成し、pickleファイル（*.pkl）に格納します。これらの系統は、sample-reconstructed strainsと呼ばれる。・・・サンプルごとに複数のsample2markers.pyスクリプトを並行して実行したい場合も、結果は同じです（クラスタシステムの設定によっては、この方法が有効な場合があります）。このステップの後、consensus_markersフォルダに全てのsample-markerファイル(*.pkl)が格納されます。

```
mkdir consensus_markers
sample2markers.py \
  -n 22 \
  --database ~/db/metaphlan/bowtie2db/221111 \
  -i sams/*.sam.bz2 \
  -o consensus_markers 
```


extract_markers.py: MetaPhlAnのデータベースから特定クレードのマーカー遺伝子配列を抽出（ex. Bacteroides_caccae (SGB1877)）。`--database`にはデータベースの「.pklファイル名」を入力すること!

-c(--clades)には基本的に種レベルゲノムビンシステム（SGB）のid（t__SGB）を指定するみたい。分類群名からSGB番号を調べるには、MetaPhlAnデータベースのspecies.txtを検索すれば良いかも。指定するSGB番号は後ろに_groupがついてる場合もある。複数指定する場合は各番号を""で囲むこと（ex. `-c "t__SGB4874" "t__SGB1836_group"`）！複数指定するとそれぞれのファイルが出力される。
- [Error when running extract_markers.py](https://forum.biobakery.org/t/error-when-running-extract-markers-py/980)
- [The bioBakery help forum - StrainPhlan4: No markers were found for the clade](https://forum.biobakery.org/t/strainphlan4-no-markers-were-found-for-the-clade/4345/3)

`-o`で指定したdb_markersフォルダにBacteroides caccae (SGB1877)のマーカー遺伝子配列が入ったファイル（「t__SGB1877.fna」など）が作成される。

```
mkdir db_markers
extract_markers.py \
  --database ~/db/metaphlan/bowtie2db/221111/mpa_vJan21_CHOCOPhlAnSGB_202103.pkl \
  -c t__SGB1877 \
  -o db_markers
```

（補足：種名からSGB番号の検索）

```
for i in s__Fusicatenibacter_saccharivorans s__Blautia_A s__Bacteroides_uniformis s__Sutterella_wadsworthensis s__Phocaeicola_dorei s__Agathobacter_faecis
do grep ${i} mpa_vJan21_CHOCOPhlAnSGB_202103_species.txt | awk '{print $1}'
done
```


アライメントと系統樹作成。output/RAxML_bestTree.t__SGB1877.StrainPhlAn4.treなどが作成される。（ファイルがけっこういっぱいできるのでcladeごとにフォルダを分けるといいかも。-oで指定するフォルダはあらかじめ作っておかないといけない。）

```
mkdir -p output/t__SGB1877
strainphlan \
  --database ~/db/metaphlan/bowtie2db/221111/mpa_vJan21_CHOCOPhlAnSGB_202103.pkl \
  -s consensus_markers/*.pkl \
  -m db_markers/t__SGB1877.fna \
  -r reference_genomes/G000273725.fna.bz2 \
  -o output/t__SGB1877 \
  -n 22 \
  -c t__SGB1877 \
  --mutation_rates
```

treeにメタデータを追加する。メタデータファイルはタブ区切りのファイルで、最初の行がヘッダー、2行目以降が各サンプルのメタデータ。メタデータは複数ファイルも指定可。output/RAxML_bestTree.t__SGB1877.StrainPhlAn4.tre.metadataが作成される（ノードラベルに|で区切ってメタデータが足される）

```
# メタデータ例（sampleIDは必須）
sampleID        subjectID
SRS055982       638754422
SRS022137       638754422
SRS019161       763496533
SRS013951       763496533
SRS014613       763840445
SRS064276       763840445
G000273725  ReferenceGenomes

# メタデータの追加
add_metadata_tree.py \
  -t output/RAxML_bestTree.t__SGB1877.StrainPhlAn4.tre \
  -f metadata.txt \
  -m subjectID \ # 使いたいメタデータヘッダー名（複数可）。デフォでは全ての列を追加する。
  --string_to_remove .fastq.bz2 # string to be removed in the tree node names
```

treやtre.metadataはRのapeなどで扱えます。