# krona

[marbl/Krona](https://github.com/marbl/Krona/wiki)

```
conda install -c bioconda krona

# セットアップ
ktUpdateTaxonomy.sh
cd /user/ngsdata/kmatsumoto/miniconda3/envs/blast/opt/krona
sh updateAccessions.sh
```

kraken2結果をかける（複数可）
```
ktImportTaxonomy \
 -q 2 \
 -t 3 \
 sample1_kraken.txt sample2_kraken.txt sample3_kraken.txt \
 -o sample1_3_kraken.html
```

blast結果（タブ）をかける
```
srun ktImportBLAST -i \
 K03851.blast K03852.blast K03852.blast K15024.blast \
 -o krona.html
```