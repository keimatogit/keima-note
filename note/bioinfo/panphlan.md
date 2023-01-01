# PanPhlAn

[SegataLab/panphlan](https://github.com/SegataLab/panphlan/wiki/Tutorial-3_0)
[macでインフォマティクス](https://kazumaxneo.hatenablog.com/entry/2017/08/07/201518)
[Distinct Genetic and Functional Traits of Human Intestinal Prevotella copri Strains Are Associated with Different Habitual Diets](https://www.sciencedirect.com/science/article/pii/S1931312819300411?via%3Dihub)

The bioBakery help forum
[Panphlan duplicate sequences](https://forum.biobakery.org/t/panphlan-duplicate-sequences/1789)

```sh
conda create -n panphlan
conda install -c bioconda panphlan

mkdir ~/db/panphlan
panphlan_download_pangenome.py -i Eubacterium_rectale -o ~/db/panphlan/
panphlan_download_pangenome.py -i Prevotella_copri -o ~/db/panphlan/
panphlan_download_pangenome.py -i Faecalibacterium_prausnitzii -o ~/db/panphlan/

conda update bowtie2

# panphlan_clean_pangenome.py を追加する
git clone https://github.com/SegataLab/panphlan.git
cp -p panphlan_clean_pangenome.py /user/ngsdata/kmatsumoto/miniconda3/envs/panphlan/bin/panphlan_clean_pangenome.py

# panphlan_clean_pangenome.pyでbiopythonが必要
conda install -c conda-forge biopython

panphlan_clean_pangenome.py --species Eubacterium_rectale --pangenome ~/db/panphlan/Eubacterium_rectale
panphlan_clean_pangenome.py --species Prevotella_copri --pangenome ~/db/panphlan/Prevotella_copri
panphlan_clean_pangenome.py --species Faecalibacterium_prausnitzii --pangenome ~/db/panphlan/Faecalibacterium_prausnitzii
```

panphlan_map.sh
```bash
#!/bin/bash

source ~/miniconda3/bin/activate panphlan

panphlan_map.py \
 -p ~/db/panphlan/${species}/${species}_pangenome.tsv \
 --indexes ~/db/panphlan/${species}/${species} \
 -i fastq/${sample}_merged.fastq.gz \
 -o map_results_${species}/${sample}_${species}.tsv
 
```
5分で終わった。

panphlan_profiling.sh
```
#!/bin/bash

source ~/miniconda3/bin/activate panphlan

panphlan_profiling.py \
 -i map_results_${species}/ \
 --o_matrix result_profile_${species}.tsv \
 -p ~/db/panphlan/${species}/${species}_pangenome.tsv \
 --add_ref

```