# awkの文字置換

6列目の「〜;g__」を「g__」に変換し、「/;s__」以降を削除。

```
head MD22074_RNA1_kraken_taxname.txt |\
awk -F "\t" '{ 
    sub( /^.+;g__/, "g__", $6)
    sub( /;s__.*$/, "", $6)
    print $2"\t"$6
}' 
```

```
# krakenのうちgenus名が付いている行のみ抽出して2,6列目を抽出、taxonomyはgenus表記だけ残して省略
cd /imetgpfs/projects/cw/IMET/221118-IMET210-mi264/02_MTT/kraken
for i in $(cat ../samples.txt)
do cat ${i}_kraken_taxname.txt | awk -F "\t" '{if(!($6 ~ /;g__;/ ) && !($6 ~ /unannotated/ )) {sub( /^.+;g__/, "g__", $6); sub( /;s__.*$/, "", $6); print $2"\t"$6} }' | sort >${i}_kraken_taxname_genus.txt
done
```