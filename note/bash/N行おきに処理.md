# N行おきに処理

## 1行めから4行おきに出力

cat file.txt | sed -ne '1~4p' 

## 3行めから4行おきに出力

cat file.txt | sed -ne '3~4p'


