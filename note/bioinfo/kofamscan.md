# kofamscan

KOGG KO番号のアノテーション。ファイルサイズが大きくてWebツールのGhostKOALAにかけるのが面倒な時に。
ひとつのクエリ配列に複数のKO番号がつくことも結構ある。その場合はひとつのクエリに対して複数行が出力される（１行あたり１KO番号）。

```sh
conda create -n kofamscan python=3.8
conda activate kofamscan
conda install -c bioconda kofamscan

>kofamscan version 1.3.0-2 has been successfully installed!
>
>This software needs a database which can be downloaded from ftp.genome.jp/pub/db/kofam/
>
>For more details see ftp://ftp.genome.jp/pub/tools/kofamscan/README
```

レファレンス：ftp://ftp.genome.jp/pub/db/kofamからprofiles.tar.gzとko_list.gzをダウンロード・解凍しておく（gunzip, tar -zxvf）。
置き場所は、~/db/GenomeNet.ftp/kofam/日付にしてある（221221）

```
cd ~/db/GenomeNet.ftp/kofam/日付
wget ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz
wget ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz # サイズ大きい
gunzip ko_list.gz
tar -zxvf profiles.tar.gz
```

ヘルプ
```
exec_annotation -h
Usage: exec_annotation [options] <query>
  <query>                    FASTA formatted query sequence file
  -o <file>                  File to output the result  [stdout]
  -p, --profile <path>       Profile HMM database
  -k, --ko-list <file>       KO information file
  --cpu <num>                Number of CPU to use  [1]
  -c, --config <file>        Config file
  --tmp-dir <dir>            Temporary directory  [./tmp]
  -E, --e-value <e_value>    Largest E-value required of the hits
  -T, --threshold-scale <scale>
                             The score thresholds will be multiplied by this value
  -f, --format <format>      Format of the output [detail]
      detail:          Detail for each hits (including hits below threshold)
      detail-tsv:      Tab separeted values for detail format
      mapper:          KEGG Mapper compatible format
      mapper-one-line: Similar to mapper, but all hit KOs are listed in one line
  --[no-]report-unannotated  Sequence name will be shown even if no KOs are assigned
                             Default is true when format=mapper or mapper-all,
                             false when format=detail
  --create-alignment         Create domain annotation files for each sequence
                             They will be located in the tmp directory
                             Incompatible with -r
  -r, --reannotate           Skip hmmsearch
                             Incompatible with --create-alignment
  --keep-tabular             Neither create tabular.txt nor delete K number files
                             By default, all K number files will be combined into
                             a tabular.txt and delete them
  --keep-output              Neither create output.txt nor delete K number files
                             By default, all K number files will be combined into
                             a output.txt and delete them
                             Must be with --create-alignment
  -h, --help                 Show this message and exit
```


## Example

インプットに圧縮ファイル(fasta.gz)は扱えなかった（21/08/19）。

同時にいくつも実行する場合（各サンプルとかで）、--tmp-dirには個別のフォルダ名を指定すること。

以下ではkofam後、出力ファイルにタブがない行の行末にタブを追加（２列のタブ区切りで読み込むため）して、ヘッダーをつけている。
--format [detail, detail-tsv, mapper, mapper-one-line]
```bash
#!/bin/bash
#SBATCH -c 96

source ~/miniconda3/bin/activate kofamscan

echo -e "\n==== kofamscan ====\n"
exec_annotation \
  gene_catalogue_aa.fasta \
  -o gene_catalogue_aa_kofam_tmp.tsv \
  --tmp-dir kofamscan_tmp \
  --cpu 96 \
  -p ~/db/GenomeNet.ftp/kofam/profiles \
  -k ~/db/GenomeNet.ftp/kofam/ko_list \
  --format mapper \
&& cat gene_catalogue_aa_kofam_tmp.tsv \
| sed -e " /\t/! s/$/\t/" \
| sed "1iContig\tkofam" > gene_catalogue_aa_kofam.tsv \
&& rm gene_catalogue_aa_kofam_tmp.tsv \
&& echo -e "\n==== END ====\n"
```
