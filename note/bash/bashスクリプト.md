# Bashスクリプト

## 名前付き変数の受け取り

get_args.sh
```
#!/bin/bash

# デフォルト値の設定
help=false
mapping=bowtie2

# 引数の受け取り
while (( $# > 0 )) # $#: 引数の個数
do
  case $1 in 
    -i | --in) # $1が"-i"か"--in"に一致する場合に以下の処理をする
        infile=$2
        shift;shift;continue;; # shift: 引数をひとつずらす, continue: ループを次に回す, ;;: 処理の終わり
    -o | --out)
        outfile=$2
        shift;shift;continue;;
    --bwa-mem )
        mapping=bwa-mem
        shift;continue;; # フラグとして使用する場合（後ろに値を取らない場合）はshift１回
    -h | --help)
        help=true
        shift;shift;break;; # break: ループを終了
    *) # 上のいずれにも当てはまらない場合
        echo "invalid option"
        exit 1
      ;;
  esac
done


# 受け取った引数を使った処理例

# ヘルプ表示
if ${help} ; then
  echo "This is help message."
  exit 1
fi

# 必須項目が空の場合の設定
if [ -z "$infile" ] || [ -z "$outfile" ]; then 
    echo "REQUIRED: --in, --out"
    exit 1
fi

# 引数の表示
echo "input file:  ${infile}"
echo "output file: ${outfile}"
echo "mapping:     ${mapping}"
```

実行例

```
./get_args.sh -i my_input -o my_output
## input file:  my_input
## output file: my_output
## mapping:     bowtie2

./get_args.sh -i my_input -o my_output --bwa-mem
## input file:  my_input
## output file: my_output
## mapping:     bwa-mem

./get_args.sh --help
## This is help message.

./get_args.sh
## REQUIRED: --in, --out

./get_args.sh -s value
## invalid option
```