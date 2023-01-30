---
---

# DADA2

DADA2はRのパッケージです。アンプリコンデータからエラーっぽい配列の除去とペアエンドのマージを行い、一塩基単位の違いに基づいたASV（amplicon sequence variant）のリスト作成してサンプル x ASVのカウント表を作成します。ASVは97% identityでのOTUなどよりも細かくグループ分けするものです。ASVには、Silvaなど参照データベースを用いて分類群アノテーションをつけることもできます（アノテーションがつかないASVもあります）。

[DADA2](https://benjjneb.github.io/dada2/)
[DADA2 Pipeline Tutorial](https://benjjneb.github.io/dada2/tutorial.html)
[filterAndTrim: Filter and trim fastq file(s).](https://rdrr.io/bioc/dada2/man/filterAndTrim.html)

[Ranking the biases: The choice of OTUs vs. ASVs in 16S rRNA amplicon data analysis has stronger effects on diversity measures than rarefaction and OTU identity threshold](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0264443)

## インストールなど

```
conda create -n dada2 python=3.8
conda install -c bioconda bioconductor-dada2
conda install -c r r-tidyverse # いっしょに使いたいので
conda install -c bioconda bioconductor-decipher bioconductor-phyloseq # いっしょに使うなら
```

確認

```
R
library(dada2)
packageVersion("dada2")
>[1] ‘1.22.0’
```

分類アノテーション用のデータセットを[Taxonomic reference data](https://benjjneb.github.io/dada2/training.html)からダウンロードしておく。
分類アノテーションでDECIFERを使う場合は[DECIFER - Downloads](http://www2.decipher.codes/Downloads.html)からダウンロードしておく（SILVA SSU r138 (modified)）。

truncLenを決めるfiraroはgithubからダウンロード（使うなら）。pythonのパッケージになっているみたい。
condaでインストールしたminiconda3/envs/dada2/bin/figaro.pyはちょっと違うみたいで、ModuleNotFoundError: No module named 'figaroSupport'というエラーが出る。


```
wget http://john-quensen.com/wp-content/uploads/2020/03/figaro.yml
conda env create -n figaro -f figaro.yml
```


```
cd
wget https://github.com/Zymo-Research/figaro/archive/master.zip
unzip master.zip
rm master.zip
mv figaro-master figaro
cd figaro/figaro
chmod 755 figaro.py
cp -p ~/figaro/figaro ~/miniconda3/envs/dada2/lib/python3.8/site-packages/ # pythonパッケージの場所に置いておく

# 確認
~/miniconda3/envs/dada2/lib/python3.8/site-packages/figaro/figaro.py --help
```

## Example

[DADA2 Pipeline Tutorial](https://benjjneb.github.io/dada2/tutorial.html)

### プライマー除去

事前にプライマーを除去しておかないといけないそうなので（プライマーのあいまい塩基のせいでキメラ配列認定されて捨てられる）、cutadaptなどで除去しておく。ついでにフィルタリングもしてしまってもいいかもと思うけど、非推奨みたい。特にリード長がバラバラになるsliding windowによる3'/5'のトリムは非推奨みたい。

```
cutadapt \
 -g AGAGTTTGATYMTGGCTCAG \
 -G GCTGCCTCCCGTAGGAGT \
 --quality-cutoff 20 \
 --minimum-length 100 \
 -o ${sample}_cutadapt_R1.fastq.gz \
 -p ${sample}_cutadapt_R2.fastq.gz \
 ${sample}_R1.fastq.gz \
 ${sample}_R2.fastq.gz
```


### 1) Getting ready （ここからR）

パッケージの呼び出しと、fastqファイルパスとサンプル名の取得

```r
library(dada2)
packageVersion("dada2") # バージョン確認
library(tidyverse) # 整形などで便利なのでtidyverseも呼び出しておく

path <- "~/mykinso/dada_test" # fastqが入っているディレクトリのパス
fnFs <- sort(list.files(path, pattern="_cutadapt_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_cutadapt_R2.fastq.gz", full.names = TRUE))
sample.names <- str_remove(basename(fnFs), "_cutadapt_R1.fastq.gz")
```

### 2) Filter and trim

リードのクオリティーフィルタリングとトリミングをして、配列データをサブディレクトリ「filtered/」に格納する。
filterAndTrimのオプションは[filterAndTrim: Filter and trim fastq file(s).](https://rdrr.io/bioc/dada2/man/filterAndTrim.html)を参照。
FigaroというツールがDADA2のtruncLenを決めるのに役立つそうです： https://github.com/Zymo-Research/figaro#figaro

```r
filtFs <- file.path(path, "filtered", paste0(sample.names, "_R1_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R2_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
out <- filterAndTrim(
    fnFs, filtFs, fnRs, filtRs,
    truncLen    = c(220,220), # 何サイクル目以降を切り捨てるか。例えば、240サイクル以降を捨てるときは、ここで240と指定します。また、これより短いreadは自動的に捨てられます。
    maxN        = 0, # Trancationを行なった後に、指定した値よりもNが多いreadを切り捨てます。DADA2はNを許さないので、無条件で０にした方がいいです。
    minQ        = 0, # Trancationを行なった後に、指定した値より低いQuality scoreを含む read を捨てます。
    maxEE       = c(2,2), # Trancationを行なった後に、指定した値より大きいexpected errorsを示した read を捨てます。expected errors：EE は sum(10^(-Q/10))によって計算されます。つまり、EE = 各readの合計のエラー率です。
    truncQ      = 2, # スコアが指定した値以下の位置以降を切り捨て？
    rm.phix     = TRUE, # 塩基の多様性を上げるために入れたphiXの配列を除きます。
    compress    = TRUE, # .gz で圧縮するかどうか。Trueで圧縮します
    multithread = FALSE # On Windows set multithread=FALSE
)
head(out)
```

filterAndTrim前にFigaroによるパラメータ決定をする場合(bash)

```bash
# ampliconLengthは「プライマー部分を除いた」ターゲット長
srun ~/miniconda3/envs/dada2/lib/python3.8/site-packages/figaro.py \
  -i ~/mykinso/01_fastq \
  -o ~/mykinso/01_fastq/figaro_res \
  --ampliconLength 310 \ 
  --forwardPrimerLength 20 \
  --reversePrimerLength 18

```


### 3) Learn the Error Rates

```r
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
```

### 4) Sample Inference

各サンプルのユニーク配列セットからASVを推定

```
dadaFs <- dada(filtFs, err=errF, multithread=TRUE) # pool=TRUEでサンプル間の情報をプール。複数のサンプルに低い頻度で存在する可能性のある配列変種に対する感度を高める。
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
# Inspecting the returned dada-class object:
dadaFs[[1]] # DADA2は〜つのユニーク配列から、〜つのバリアントを推定しました。
```

### 5) Merge paired reads

デフォルトでは、オーバーラップが12塩基以上で100% identity（引数で変更可）。
ほとんどのリードが正常にマージされるはずで、そうでない場合は上流過程のパラメータを見直す必要があるかもしれません。filterAndTrimのtruncLenでオーバーラップ部分を削除してしまったとか？

```r
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
```

### 6) Construct sequence table

行がサンプル、列がASVのカウント表を作成する。ASVで予想よりはるかに長い配列や短い配列はああった場合は、非特異的なプライミングの結果である可能性があるので除く。

```r
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# ASVの配列長の分布
table(nchar(getSequences(seqtab)))

# 以下でターゲット以外の長さの配列を配列テーブルから削除できます。
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 280:350]
dim(seqtab2)
```

### 7) Remove chimerasRemove chimeras

キメラ配列の頻度はデータセットによって大きく異なり、実験手順やサンプルの複雑さなどの要因に依存します。Tutrialの例では、全ASVのうち21%、存在量を考慮した場合は（merged 配列のうちのどのくらいかは）4%。基本的にはキメラ除去後もほとんどのリードが残っているはずです（ただし、配列バリアントの大部分が除去されることも珍しくありません）。キメラとして除去されたリードが多い場合、上流過程の見直しが必要な場合があります。ほとんどの場合、DADA2パイプラインを開始する前に除去されなかった**曖昧な塩基を持つプライマー配列**が原因です。

```r
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# 全ASVのうちキメラ配列の割合
(ncol(seqtab) - ncol(seqtab.nochim))/ncol(seqtab)

# 存在量を考慮したキメラ配列の割合（merged配列のうちどのくらいか）
(sum(seqtab) - sum(seqtab.nochim))/sum(seqtab)
```

### 8) Track reads through the pipeline

各ステップを通過したリードの数をチェック！
- フィルタリング以外では、大部分のリードが失われるようなステップはないはずである。
- 大部分のリードがマージに失敗した場合、フィルタリングのステップで使用したtruncLenパラメータを再検討し、切り捨てられたリードがアンプリコンにまたがることを確認する必要があるかもしれません。
- 大部分のリードがキメラとして除去された場合、プライマーの除去を再検討する必要があるかもしれません。除去されていないプライマーのあいまいなヌクレオチドがキメラの同定に干渉するためです。

```r
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track
```

### 9) Assign taxonomy

ナイーブベイズ分類。

```r
# 参照DBと逆向きになっている場合は、`tryRC=TRUE`
taxa <- assignTaxonomy(seqtab.nochim, "~/db/dada2/silva138.1/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE) # 97% identity?
taxa <- addSpecies(taxa, "~/db/dada2/silva138.1/silva_species_assignment_v138.1.fa.gz") # 100% identity?

# Let’s inspect the taxonomic assignments:
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# アノテーションがついたASV数の確認
nrow(taxa.print) # 全ASV数
sum(!is.na(as.data.frame(taxa.print)$Genus))   # GenusがついたASV数
sum(!is.na(as.data.frame(taxa.print)$Species)) # SeciesがついたASV数
```

### 9) Assign taxonomy (DECIPHER)

IDTAXAアルゴリズムというのを使っているらしい。assignTaxonomyの代替案。

```r
dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
load("~/db/decifer/SILVA_SSU_r138_2019.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest

# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
        m <- match(ranks, x$rank)
        taxa <- x$taxon[m]
        taxa[startsWith(taxa, "unclassified_")] <- NA
        taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)

# Let’s inspect the taxonomic assignments:
taxid.print <- taxid # Removing sequence rownames for display only
rownames(taxid.print) <- NULL
head(taxid.print)

# アノテーションがついたASV数の確認
nrow(taxid.print) # 全ASV数
sum(!is.na(as.data.frame(taxid.print)$genus))   # GenusがついたASV数
sum(!is.na(as.data.frame(taxid.print)$species)) # SeciesがついたASV数
```

### Appendix) Write tables

テーブルの出力（ASVには連番の名前をつける）。taxaのrownamesとseqtab.nochimのcolnamesが同じ（ASVの配列）なことに注意。（taxaはASV x 分類ランク、seqtab.nochimはサンプル x ASVのカウント）

```
asv_names <- paste0("ASV", str_pad(1:nrow(taxa), nchar(nrow(taxa)), pad=0))
asv_seq  <- tibble(asv = asv_names, seq = rownames(taxa))
asv_taxa <- tibble(asv = asv_names) %>%
    bind_cols(as_tibble(taxa)) %>%
    mutate(taxonomy = paste0("k__", Kingdom, ";p__", Phylum, ";c__", Class, ";o__", Order, ";f__", Family, ";g__", Genus, ";s__", Species))
asv_count <- t(seqtab.nochim) %>%
    as_tibble() %>%
    mutate(asv = asv_names) %>%
    select(asv, everything())
write_tsv(asv_seq,   paste0(path, "/ASV_seq.tsv")) # ASVの配列
write_tsv(asv_taxa,  paste0(path, "/ASV_taxonomy.tsv")) # ASVの分類アノテーション
write_tsv(asv_count, paste0(path, "/ASV_count.tsv")) # ASV x サンプルのカウントデータ
```

### Appendix) Handoff to phyloseq

DADA2の結果をphyloseqに渡してプロットなど出す場合の例

ライブラリ読み込みなど

```r
library(phyloseq)
library(Biostrings)
library(ggplot2)
theme_set(theme_bw())
```

サンプルメタデータの作成

```r
samples.out <- rownames(seqtab.nochim)
person <- samples.out %>%
    str_split("_", simplify = TRUE) %>%
    .[,1]
day <- samples.out %>%
    str_split("_", simplify = TRUE) %>%
    .[,2] %>%
    as.Date("%y%m%d")
samdf <- data.frame(Person = person, Day = day)
rownames(samdf) <- samples.out
```

dada2の出力からphyloseqオブジェクトを構築。phyloseqではASVの完全なDNA配列よりも短い名前（例えばASV21）を使用したいけど、DNA配列もとっておきたい。なので、ASVのDNA配列をphyloseqオブジェクトのrefseqスロットに格納し、ASV名を短い文字列に変更している。

```r
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
colnames(ps@otu_table) # ASV[連番]になってる
```

Visualize alpha-diversity:

```r
pdf("~/mykinso/plot_richness.pdf")
plot_richness(ps, x="Day", measures=c("Shannon", "Simpson"), color="Person")
dev.off()
```

Ordinate (beta-diversity):

```r
# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

pdf("~/mykinso/plot_ordinate.pdf")
plot_ordination(ps.prop, ord.nmds.bray, color="Person", title="Bray NMDS")
dev.off()
```

Bar plot:

```r
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
pdf("~/mykinso/plot_top20.pdf")
plot_bar(ps.top20, x="Day", fill="Family") + facet_wrap(~Person, scales="free_x") + scale_fill_brewer(palette="Set2")
dev.off()
```