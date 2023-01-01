# コマンドラインでのプログラム実行時に引数を扱う

<!-- @import "[TOC]" {cmd="toc" depthFrom=2 depthTo=6 orderedList=false} -->
<!-- code_chunk_output -->
v
- [python](#python)
- [R](#r)
- [shell](#shell)

<!-- /code_chunk_output -->


## python

argparseライブラリを使う。

[argparse](https://docs.python.org/3/library/argparse.html)
[ArgumentParserの使い方を簡単にまとめた](https://qiita.com/kzkadc/items/e4fc7bc9c003de1eb6d0)


```python

import argparse

parser = argparse.ArgumentParser(
    description     = 'Description of this program (help message).',
    formatter_class = argparse.ArgumentDefaultsHelpFormatter,
    add_help        = False # 省略形'-h'を使用したい場合、-hによるヘルプ表示を無効にしておく
)

# Positional arguments（位置指定引数。並べた順番に引数が引き渡される。入力必須になる。）
parser.add_argument('name', help='name')
parser.add_argument('num',  help='number', type=int)

# Optional arguments（名前つき引数。required=Trueで入力必須になる。）
parser.add_argument('-i', '--infile',  help='input  file name. Required.', required=True)
parser.add_argument('-o', '--outfile', help='output file name. Required.', required=True)
parser.add_argument('-h', '--host',    help='host name', default='localhost')

# 上でヘルプ表示を無効にした場合も--helpで表示できるようにする
parser.add_argument('--help', action='help')


args = parser.parse_args()

print('name:    ' + args.name)
print('number:  ' + str(args.num))
print('infile:  ' + args.infile)
print('outfile: ' + args.outfile)
print('host:    ' + args.host)
```

実行例
```
python3 commandargs.py myName 2 -i infile.txt -o outfile.txt
```

結果
```
name:    myName
number:  2
infile:  infile.txt
outfile: outfile.txt
host:    hostname
```

引数が多い場合、名前つき引数の方が順番を覚えなくて良いので良いと思う。


## R

argparseパッケージを使う。使い方はpythonのargparseライブラリと基本的に同じ。

[trevorld/r-argparse - GitHub](https://github.com/trevorld/r-argparse)

```r

library(argparse)

parser <- ArgumentParser(
  description = "Load csv data into mariadb",
  add_help    = FALSE # 省略形'-h'を使用したい場合、-hによるヘルプ表示を無効にしておく
)

# Positional arguments（位置指定引数。並べた順番に引数が引き渡される。入力必須になる。）
parser$add_argument("name", help="name")
parser$add_argument("num",  help="number", type="integer")

# Optional arguments（名前つき引数。required=TRUEで入力必須になる。）
parser$add_argument("-i", "--infile",  help="input  file name. Required.", required=TRUE)
parser$add_argument("-o", "--outfile", help="output file name. Required.", required=TRUE)
parser$add_argument("-h", "--host",    help="host name", default="localhost")

# ヘルプ表示を追加
parser$add_argument("--help", action="help")


args <- parser$parse_args()

cat("name:    ", args$name,    "\n")
cat("number:  ", args$num,     "\n")
cat("infile:  ", args$infile,  "\n")
cat("outfile: ", args$outfile, "\n")
cat("host:    ", args$host,    "\n")
```

実行例
```
Rscript commandargs.R myName 2 -i infile.txt -o outfile.txt
```

結果
```
name:     myName
number:   2
infile:   infile.txt
outfile:  outfile.txt
host:     localhost
```

## shell

位置指定引数のみ。名前付き引数の方法は分からない。

```shell:commandargs.sh

#!/bin/sh
#$ -S /bin/sh

name=$1 #'='の前後にスペースを入れてはいけない（変数名ではなく文字列として認識される）
num=$2
infile=$3
outfile=$4
host=${5:-localhost}

echo ${name}
echo ${num}
echo ${infile}
echo ${outfile}
echo ${host}
```

実行例
```
sh commandargs.sh myName 2 infile.txt outfile.txt
```

結果
```
myName
2
infile.txt
outfile.txt
localhost
```
