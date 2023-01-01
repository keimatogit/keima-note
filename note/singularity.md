# singularity


## ubuntuにsingularityをインストール
こちらに沿って行います
[Quick Start](https://sylabs.io/guides/3.4/user-guide/quick_start.html)
WSL1だとうまくいかなかったので、WSL2にしておきましょう(2021/05/24)


1. 必要なライブラリをインストール
```
sudo apt-get update && sudo apt-get install -y \
  build-essential \
  libssl-dev \
  uuid-dev \
  libgpgme11-dev \
  squashfs-tools \
  libseccomp-dev \
  wget \
  pkg-config \
  git \
  cryptsetup
```

2. goをインストール
[GO install](https://golang.org/doc/install)
古いものがある場合は先にrmしておく。
```
wget https://golang.org/dl/go1.16.4.linux-amd64.tar.gz
[[sudo] rm -rf /usr/local/go &&] [sudo] tar -C /usr/local -xzf go1.16.4.linux-amd64.tar.gz
# export PATH=$PATH:/usr/local/go/bin を.basrcや.profileに追加
go version # 確認
```

3. singularityのインストール
[singularity releases](https://github.com/sylabs/singularity/releases)
```
wget https://github.com/sylabs/singularity/releases/download/v3.7.3/singularity-3.7.3.tar.gz
tar -xzf singularity-3.7.3.tar.gz
cd singularity
./mconfig && \
  make -C builddir && \
  sudo make -C builddir install
singularity --version # 確認
```

## コンテナイメージのダウンロード

Singularity Library からダウンロード
```
singularity pull library://path
```

Docker Hub からダウンロード
```
singularity pull docker://path
```
例
```
singularity pull docker://jupyter/datascience-notebook
```

## コンテナを使う

コンテナに入ってその中のシェルを使う
```
singularity shell my_container.sif
```
直接コマンドを実行させる
```
singularity exec my_container.sif command
```
例
```
singularity exec datascience-notebook_latest.sif python
```
## マウント

[Bind Paths and Mounts](https://sylabs.io/guides/3.0/user-guide/bind_paths_and_mounts.html)

以下の場所は自動でマウントされ、コンテナ内からアクセスできます。
- $HOME
- $PWD
- /tmp
- /proc
- /sys
- /dev

他の場所にアクセスしたい場合は、実行時に `--bind`オプションで指定するか、シェルの環境変数SINGULARITY_BINDPATHを定義すると、`/mnt`にマウントされます。
```
# ホストが/optをコンテナの/optに、ホストの/dataがコンテナの/mntにマウントされる
singularity shell --bind /opt,/data:/mnt my_container.sif
```
```
# 環境変数を定義（いつも使う場合、.bashrcなどに書き込んでおく）
export SINGULARITY_BIND="/opt,/data:/mnt"
```
