# Prodigal

[Wiki Documentation](https://github.com/hyattpd/prodigal/wiki)

```
conda install -c bioconda prodigal
```

```
prodigal -v

Prodigal V2.6.3: February, 2016
```


## ヘルプ

```
>prodigal -h

Usage:  prodigal [-a trans_file] [-c] [-d nuc_file] [-f output_type]
                 [-g tr_table] [-h] [-i input_file] [-m] [-n] [-o output_file]
                 [-p mode] [-q] [-s start_file] [-t training_file] [-v]

         -a:  Write protein translations to the selected file.
         -c:  Closed ends.  Do not allow genes to run off edges.
         -d:  Write nucleotide sequences of genes to the selected file.
         -f:  Select output format (gbk, gff, or sco).  Default is gbk.
         -g:  Specify a translation table to use (default 11).
         -h:  Print help menu and exit.
         -i:  Specify FASTA/Genbank input file (default reads from stdin).
         -m:  Treat runs of N as masked sequence; don't build genes across them.
         -n:  Bypass Shine-Dalgarno trainer and force a full motif scan.
         -o:  Specify output file (default writes to stdout).
         -p:  Select procedure (single or meta).  Default is single.
         -q:  Run quietly (suppress normal stderr output).
         -s:  Write all potential genes (with scores) to the selected file.
         -t:  Write a training file (if none exists); otherwise, read and use
              the specified training file.
         -v:  Print version number and exit.
```


## 使い方

オプション`-i`または標準入力で配列を渡す。圧縮ファイルの場合は`gunzip -c`で渡してやるといい。
```
gunzip -c sample.fasta.gz |
prodigal \
  -p meta \
  -g 11 \
  -o sample.gbk \
  -d sample.orf.fasta \
  -a sample.aa.fasta
```

長さは20aaから出力されている（V2.6.3）。min lenの設定はないみたい？

出力ファイルの配列は折り返されているので、配列中の改行を除きたい場合は`seqkit seq -w 0 sample.aa.fasta >sample.aa.w0.fasta`など。

各サンプルからのorf配列名にサンプル名を入れてcd-hitにかける用にまとめる場合
```
for i in $(cat samples.txt)
do seqkit seq -w 0 ${sample}.orf.fasta | sed -e "1~2s/ #.*$//" | \
  sed -e "1~2s/^>/>${sample}_/" >> all_orf.fasta
```

## 出力

### orfやアミノ酸配列のfasta

[Understanding the Prodigal Output](https://github.com/hyattpd/prodigal/wiki/understanding-the-prodigal-output)
- 元のリード名_連番（これがユニークであるかは元のリード名による）
- # 開始位置
- # 終了位置
- # ストランド（1か-1）
- # ID（リード連番_リード内の連番のユニークなID）
- ;partial: 遺伝子が完全かどうか、0が境界ありで1が境界なし（未完成）。00:開始と終止コドンをもつ、01:開始コドンのみ、10:終止コドンのみ、11:どちらもなし。
- ;start_type: 開始コドン（開始コドンがない場合"Edge"）
- ;stop_type:  停止コドン（停止コドンがない場合"Edge"）
- ;rbs_motif: Prodigalが見つけたRBS(ribosomal binding site)モチーフ
- ;rbs_spacer: 開始コドンとRBSモチーフの間の長さ
- ;gc_cont: GC content



## prodigalの翻訳結果とtranseqやseqkit tranlateの翻訳結果との違い

prodigalで冒頭TTGがMだったところが、seqkit_translateやseqkit translateでLになっていた（冒頭TTGが全部そうかは確認してない）。genetic codeはどちらも11を指定。（seqkit translateの --init-codon-as-MでMになる）

```sh
grep "k141_2_1" -3 prodigal.orf.fasta
>k141_2_1 # 2 # 307 # -1 # ID=16_1;partial=10;start_type=TTG;rbs_motif=GGA/GAG/AGG;rbs_spacer=5-10bp;gc_cont=0.415
TTGGCAGTTATACAGACAGGAAGGGGTGTACCGTATATGAGACAGATCATAGTAGAAGGACAGACTCTGG

(bio) bash-4.2$ sed -n 60p prodigal.aa.fasta
MAVIQTGRGVPYMRQIIVEGQTLDVKDGTTYLELAKNFQKKFDHDIVLVLENNKMRELFR
(bio) bash-4.2$ sed -n 60p seqkit_translate.fasta
LAVIQTGRGVPYMRQIIVEGQTLDVKDGTTYLELAKNFQKKFDHDIVLVLENNKMRELFR
grep "k141_2_1" -3 transeq.fasta
LAVIQTGRGVPYMRQIIVEGQTLDVKDGTTYLELAKNFQKKFDHDIVLVLENNKMRELFR
```
