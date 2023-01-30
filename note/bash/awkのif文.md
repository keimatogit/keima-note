# awk内でbash変数を使う

```
# 14列目が指定のパターンを含む場合12列目をプリント
cat filereport_PRJEB38984.txt | sed '1d' | awk -F "\t" '{ if($14 ~ /^W000[123456789]_/ || $14 ~ /^W0011/ ) {print $12} }' 
```


6列目が「^sk__Bacteria;」を含み、かつ「;g__;」を含まない場合に2列目と6列目をプリント
```
awk -F "\t" '{ if($6 ~ /^sk__Bacteria;/ && !($6 ~ /;g__;/ )) {print $2"\t"$6} }'
```

```
# kofamのうちマーカーKOの行だけを抽出（MD22083_RNA2はまだ終わってない）
cd /imetgpfs/projects/cw/IMET/221118-IMET210-mi264/02_MTT/bbmerge/kofam
for i in MD22074_RNA1 MD22074_RNA2 MD22083_RNA1
do cat ${i}_kofam.tsv | awk -F "\t" '{if($2!="") {print $0} }' | awk -F "\t" '{if( $2=="K06942" || $2=="K01889" || $2=="K01887" || $2=="K01875" || $2=="K01883" || $2=="K01869" || $2=="K01873" || $2=="K01409" || $2=="K03106" || $2=="K03110" ) {print $0} }' | sort >${i}_kofam_markerKOs.tsv
done

# krakenのうちgenus名が付いている行のみ抽出して2,6列目を抽出、taxonomyはgenus表記だけ残して省略
cd /imetgpfs/projects/cw/IMET/221118-IMET210-mi264/02_MTT/kraken
for i in $(cat ../samples.txt)
do cat ${i}_kraken_taxname.txt | awk -F "\t" '{if(!($6 ~ /;g__;/ ) && !($6 ~ /unannotated/ )) {sub( /^.+;g__/, "g__", $6); sub( /;s__.*$/, "", $6); print $2"\t"$6} }' | sort >${i}_kraken_taxname_genus.txt
done
```