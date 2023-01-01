# InStrain

メタゲノムのSNVs(一塩基変異)の探索。入力はレファレンスゲノムデータベースにマッピングしたsamまたはbamファイル。


[github](https://github.com/MrOlm/inStrain)
[InStrain](https://instrain.readthedocs.io/en/latest/)

[macでインフォマティクス - メタゲノムデータから集団の微細多様性をプロファイリングする inStrain](https://kazumaxneo.hatenablog.com/entry/2022/05/05/204219)

## インストール

v1.6.3（2022-11-31）

```
conda create -n instrain
conda activate instrain
conda install -c conda-forge -c bioconda -c defaults instrain

# 確認
inStrain --help
```

## レファレンスゲノムの準備

[inStrain - Tutrial#2](https://instrain.readthedocs.io/en/latest/tutorial.html#tutorial-2-running-instrain-using-a-public-genome-database)
[Important concepts - 4. Establishing and evaluating genome databases](https://instrain.readthedocs.io/en/latest/important_concepts.html#establishing-and-evaluating-genome-databases)


### レファレンスについて

必要に応じて、
- 既存データベースを使う
- アセンブル（カスタム）ゲノムを使う
- アセンブル（カスタム）ゲノム + 既存データベースを使う
ことを考えましょう。

マッピングを行う際には、すべてのゲノムに同時にマッピングすることが重要です。マッピングしたい全てのゲノムを1つの.fastaファイルにまとめましょう。

メタゲノムのマップ率は、ゲノムデータベースの充実具合の目安となる。土壌では20-40％、ヒトマイクロバイオームでは60-80％、単純でよく定義されたコミュニティでは90％以上と予想されるとのこと。


### 既存データベースからレファレンス準備（human gut prokaryotes)


[inStrain - Tutorial #2) Running inStrain using a public genome database](https://instrain.readthedocs.io/en/latest/tutorial.html#tutorial-2-running-instrain-using-a-public-genome-database)

[mgnify_genomes](http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/)にはhuman-gutやhuman-oralやmarineのゲノムデータベースがある。FTPサイトは[こちら](http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/)。新しいバージョンからダウンロードしよう。

UHGG (Unified Human Gastrointestinal Genome)はヒト腸内原核生物の非冗長ゲノムセット（[Almeida et al. 2021](https://www.nature.com/articles/s41587-020-0603-3)）。これをダウンロードする。v2.0.1のREADMEによると、4,744種にクラスタリングされる289,232ゲノムが含まれてるらしい。

まずはメタデータをダウンロード

```
cd ~/db/mgnify_genomes/human-gut/v2.0.1 # ダウンロードしたい場所に移動
wget http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.1/genomes-all_metadata.tsv
```

メタデータをもとに１種につき１ゲノム配列をダウンロード
inStrainのTutrialを参照して、v2.0.1のメタデータの形式に合わせてbash用に変え、シェルスクリプトにしてsbatchしています。

```
cat genomes-all_metadata.tsv | awk -F "\t" '{if ($14 == $1) print $14}' >species_rep_list.txt
wc -l species_rep_list.txt # 行数が種数と同じことを確認

# シェルスクリプトを作成
for sp_accession in $(cat species_rep_list.txt)
do echo "wget http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.1/species_catalogue/${sp_accession:0:11}/${sp_accession}/genome/${sp_accession}.fna" >> UHGG_reps/wget.sh
done

# wget.shの冒頭にシバン等を追加して、sbatchで実行
cd UHGG_reps
sbatch -J get_UHGG_reps -o slurm_wget_%j.out ./wget.sh

# 種数と同じ数のfastaファイルがダウンロードできてるか確認
ls UHGG_reps/*.fna | wc -l
```

scaffold-to-bin fileというものを作成する。２列のタブ区切りテキストで、1列目はscaffoldの名前、2列目はscaffoldが属するゲノムファイル名。

```
sbatch -p c6420 -c 22 -J parse_stb -o slurm_parse_stb.out --wrap="parse_stb.py --reverse -f UHGG_reps/*.fna -o UHGG_reps.stb"
```

Prodigalで各ゲノムから遺伝子予測をして、各ゲノムの予測遺伝子ファイルを作成する。

```
mkdir UHGG_genes
cd UHGG_reps/

# シェルスクリプトを作成
for sp_accession in $(ls *.fna | sed -e "s/.fna$//g")
do echo "prodigal -p single -m -i ${sp_accession}.fna -o ../UHGG_genes/${sp_accession}.genes -a ../UHGG_genes/${sp_accession}.gene.faa -d ../UHGG_genes/${sp_accession}.gene.fna" >>prodigal.sh
done 

# prodigal.shにシバンや「source ~/miniconda3/bin/activate instrain」を書き加えて、sbatchにかける
sbatch -J UHGG_prodigal -o slurm_prodigal_%j.out ./prodigal.sh

# 種数分出力されていることを確認
ls UHGG_genes/*fna | wc -l
ls UHGG_genes/*faa | wc -l
ls UHGG_genes/*genes | wc -l
```

各ゲノムの予測遺伝子を１ファイルにまとめる。

```
cat UHGG_genes/*.gene.fna > UHGG_reps.genes.fna
cat UHGG_genes/*.gene.faa > UHGG_reps.genes.faa
```

各ゲノムのfastaファイルも1ファイルにまとめておく。

```
cat UHGG_reps/*.fna > UHGG_reps.fasta
```


## マッピング

マッピングインデックスの作成（マッピングソフトは任意だけど今回はbowtie2）
マニュアルでは`--large-index`がついてたので付けておく
(221207 bowtie2 v2.4.4で作成)

```
conda activate bio
sbatch -p c6420 -c 22 -J bowtie2-build -o slurm_bowtie2-build_%j.out --wrap="\
bowtie2-build UHGG_reps.fasta ~/db/bowtie2/UHGG_reps/UHGG_reps --large-index --threads 22"
```

マッピング
マップされなかったリードはinStrainでは使われない。マップ率を確認しておこう。

```
bowtie2 \
 -x ~/db/bowtie2/UHGG_reps/UHGG_reps \
 -1 ../01_bowtieToHuman/${sample}_unmap_R1.fastq.gz \
 -2 ../01_bowtieToHuman/${sample}_unmap_R2.fastq.gz \
 -p 22 \
 -S ${sample}.sam
```

## Running inStrain profile

> When running inStrain on a big database like we have here it’s critical to add the flag --database mode. This flag does some quick calculations to figure out which genomes are probably not present, and stops working on them right away. This leads to dramatic reductions in RAM usage and computational time.

```
inStrain profile \
  ${sample}.sam \
  ~/db/mgnify_genomes/human-gut/v2.0.1/UHGG_reps.fasta \
  -o ${sample}.IS \
  -p 10 \
  -g ~/db/mgnify_genomes/human-gut/v2.0.1/UHGG_reps.genes.fna \
  -s ~/db/mgnify_genomes/human-gut/v2.0.1/UHGG_reps.stb \
  --database_mode
```
${sample}.ISフォルダ内に４つのフォルダができる。
figures  log  output  raw_data

## 'inStrain profile'の出力

[Expected output](https://instrain.readthedocs.io/en/latest/example_output.html)を参照。

outputフォルダ内
- *gene_info.tsv
- *linkage.tsv
- *scaffold_info.tsv
- *genome_info.tsv
- *mapping_info.tsv
- *SNVs.tsv.gz

### （例）output/*scaffold_info.tsv

サンプルで見つかったレファレンスscaffoldごとの基本的な情報

- scaffold: レファレンスのscaffold配列名
- scaffold: レファレンスのscaffold配列の長さ
- coverage: scaffoldの各塩基（位置）のカバレッジの平均値
- breadth: 少なくとも１リードでカバーされているscaffold内の塩基の割合（1だと全領域がカバーされている）
- nucl_diversity: 塩基の多様性があった箇所の[nucleotide diversity](https://instrain.readthedocs.io/en/latest/overview.html#term-nucleotide-diversity)の平均値。つまり塩基の多様性があった位置が一箇所のみの場合は、その箇所のnucleotide diversityの値となる。min_covを満たす箇所がない場合は空白となる。
- coverage_median: カバレッジの中央値（カバレッジ0の位置も含めた時の中央値）
- coverage_std: カバレッジの標準偏差
- coverage_SEM: カバレッジの平均値の標準誤差
- breadth_minCov: min_cov以上のカバレッジを持つ塩基の割合
- breadth_expected: breadsの期待値（scaffold全体が均等に含まれる場合の期待値をカバレッジの値から得ている）。実際のbreadthの値がこれよりすごく小さい場合は、トランスポゾンやプロファージなどの特定の領域にのみマッピングされていることを示す。
- nucl_diversity_median: 塩基多様性が計算できた位置の塩基多様性の中央値
- nucl_diversity_rarefied: rarefied_coverage (デフォルトでは50) 以上のカバレッジがある位置の平均塩基多様性
- nucl_diversity_rarefied_median: 上と同様に中央値
- breadth_rarefied: rarefied_coverage以上のカバレッジがある位置のパーセンテージ
- conANI_reference: Consensus ANI。レファレンス配列とサンプルリードの平均塩基同一性。通常のANIに相当。マイナーアレルは無視される。
- popANI_reference: Population ANI。inStrain独自のANIでマイナーアレルを無視しない。例えばレファレンス配列がAの位置でマップされたリードの60本がCで40本がAの場合、違う箇所としてカウントされない（conANIではカウントされる）。
- SNS_count: Single nucleotide substitution（一塩基置換）の個数。レファレンス配列でAの位置で、すべてのリードはその位置がCだった場合は、SNSにカウントされる（リードの半分がAでもう半分がCの場合は、SNVにカウントされる）。
- SNV_count: Single nucleotide variant（一塩基変異）の個数。リードセット内で遺伝的変異のある位置の個数を示す。
- divergent_site_count: SVSまたはSNVがあった箇所の個数。
- consensus_divergent_sites:
- population_divergent_sites: 


## Running inStrain compare

inStrain profileの結果をサンプル間で比較する。

```
inStrain compare \
  -i ${sample1}.IS/ \
     ${sample2}.IS/ \
  -s ~/db/mgnify_genomes/human-gut/v2.0.1/UHGG_reps.stb \
  -p 6 \
  -o UHGG_reps.IS.COMPARE \
  --database_mode
```

出力ファイルについては、[Expected output - inStrain compare](https://instrain.readthedocs.io/en/latest/example_output.html#instrain-compare)を参照