# Figaro 

QIIMEやDADA2でのtruncLenなどのパラメータ決定。まだ使えていない。

[John Quensen - FIGARO](https://john-quensen.com/tutorials/figaro/)

## インストール

普通に`conda install -c bioconda figaro`では使えなかった（2022/10/14）。

[John Quensen - FIGARO](https://john-quensen.com/tutorials/figaro/)

```
wget http://john-quensen.com/wp-content/uploads/2020/03/figaro.yml
conda env create -n figaro -f figaro.yml

cd
wget https://github.com/Zymo-Research/figaro/archive/master.zip
unzip master.zip
rm master.zip
mv figaro-master figaro
cd figaro/figaro
chmod 755 figaro.py

cp -pr ~/figaro/figaro ~/miniconda3/envs/figaro/lib/python3.6/ # pythonパッケージの場所に置いておく

~/miniconda3/envs/figaro/lib/python3.6/figaro/figaro.py --help
```

↑でヘルプまでは出たけど、実際のデータで実行しようとすると色々エラーが出て使えなかった（2022/10/14, fastqHandler.FastqValidationError: Unable to validate fastq files enough to perform this operation. Please check log for specific error(s).）。


```
# ampliconLengthは「プライマー部分を除いた」ターゲット長
~/miniconda3/envs/figaro/lib/python3.6/figaro/figaro.py \
  -i ~/mykinso/01_fastq \
  -o ~/mykinso/01_fastq/figaro_res \
  --ampliconLength 310 \
  --forwardPrimerLength 20 \
  --reversePrimerLength 18
```