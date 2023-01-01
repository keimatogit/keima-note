# transeq（EMBOSS: transeq）

ORF予測。

[emboss transeq](https://www.bioinformatics.nl/cgi-bin/emboss/help/transeq)

インストール
```
conda install -c bioconda emboss
```

## 使い方

- 配列名中に"|"が入っていると名前を変えられてしまうよ！
- 圧縮ファイル（.gz）は扱えない

```sh
gunzip -c my_orf_catalogue.fasta.gz > my_orf_catalogue.fasta

srun transeq \
  -sequence my_orf_catalogue.fasta \
  -outseq   my_orf_catalogue_aa_tmp.fasta \
  -frame    1 \
  -table    11

# 配列中の改行と配列名の最後についた_1(frame)を外しておく
srun cat my_orf_catalogue_aa_tmp.fasta | \
seqkit replace -p "_1$" | \
seqkit seq -w 0 > my_orf_catalogue_aa.fasta

gzip my_orf_catalogue_aa.fasta
rm my_orf_catalogue.fasta my_orf_catalogue_aa_tmp.fasta
```

.gzが扱えないので、`seqkit translate`でもいいかも。
参考：`seqkit translate --frame 1 --transl-table 11 input.fasta.gz > output.fa`


