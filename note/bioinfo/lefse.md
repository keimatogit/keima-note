# LEfSe


## Install

[LEfSe](https://huttenhower.sph.harvard.edu/lefse/)
[LEfSeの使用法](https://note.com/ytomy/n/n39b64afd3862)

[condaでLEfSeを導入するときに注意すること](https://qiita.com/FlaskSakaguchi/items/b9e70129793174d939e1)

```
conda create -n lefse python=3.7 -y
conda activate lefse
conda install -c bioconda -c conda-forge lefse -y
```


`conda create -n lefse -c bioconda lefse`だと`lefse_run.py`でエラーが出た（pythonのバージョンの問題？）

```
lefse_run.py
# Traceback (most recent call last):
#   File "/user/ngsdata/kmatsumoto/miniconda3/envs/lefse/lib/python3.10/site-packages/rpy2/rinterface_lib/openrlib.py", line 17, in <module>
#     import _rinterface_cffi_api as _rinterface_cffi  # type: ignore
# ImportError: libffi.so.7: cannot open shared object file: No such file or directory
```


## Input format

> The text tab-delimited input file consists of a list of numerical features, the class vector and optionally the subclass and subject vectors. The features can be read counts directly or abundance floating-point values more generally, and the first field is the name of the feature. Class, subclass and subject vectors have a name (the first field) and a list of non-numerical strings.

> Although both column and row feature organization is accepted, given the high-dimensional nature of metagenomic data, the listing of the features in rows is preferred.

- タブ区切りの数値表。
- １列目はラベル（文字列）。
- １行目がグループ名、２行目？がサブグループ名（任意）、３行目がID（これらも１列目に名称をつける）。４行目以降が特徴量。
- 数値はカウントでも頻度でも。
- 行が分類群や遺伝子で、列がサンプル。
- 分類群名は"|"で区切ると上の分類群も検定してくれる。（上の分類群の行がなくても勝手に集計してくれるもよう。行を入れる場合は合計頻度を入れるもよう（サンプルファイル参照）。）。"|"で区切るとcladogramも出せる。これらをしてほしくない時は|で区切らなければOK。
- 比較は１行目のグループ間で行われる。サブグループは、サブグループごとでも（それぞれのサブグループ内でも）グループ間で違うか、というように使われるみたい。[LefSe Class and Subclass](https://groups.google.com/g/lefse-users/c/QPT7S-o-RhM)


例）[sample input file: hmp_small_aerobiosis.txt](https://github.com/biobakery/biobakery/raw/master/test_suite/biobakery_tests/data/lefse/input/hmp_small_aerobiosis.txt)
```
oxygen_availability	High_O2	Mid_O2	Low_O2	Mid_O2	Low_O2
body_site	ear	oral	gut	oral	gut
subject_id	158721788	158721788	159146620	159005010	159166850
Archaea|Euryarchaeota|Methanobacteria|Methanobacteriales|Methanobacteriaceae|Methanobrevibacter	2.96541e-06	5.08937e-06	4.93921e-06	7.86782e-06	6.22604e-06
Bacteria	0.999994	0.99999	0.99999	0.999984	0.999988
Bacteria|Acidobacteria	5.0412e-05	8.65194e-05	8.39666e-05	0.000133753	0.000105843
Bacteria|Acidobacteria|Acidobacteria_Gp10|Gp10	2.96541e-06	5.08937e-06	4.93921e-06	7.86782e-06	6.22604e-06
Bacteria|Acidobacteria|Acidobacteria_Gp11|Gp11	2.96541e-06	5.08937e-06	4.93921e-06	7.86782e-06	6.22604e-06
Bacteria|Acidobacteria|Acidobacteria_Gp16|Gp16	2.96541e-06	5.08937e-06	4.93921e-06	7.86782e-06	6.22604e-06
Bacteria|Acidobacteria|Acidobacteria_Gp17|Gp17	2.96541e-06	5.08937e-06	4.93921e-06	7.86782e-06	6.22604e-06
```


## ヘルプメッセージ

### lefse_format_input.py -h

- normalizationは、サンプル（列）ごとに列の合計値で割った後に`-o`の値を掛ける？[Questions about pre-sample normalization and low aboundant taxa for LEfSe](https://groups.google.com/g/lefse-users/c/5IFgZRjHGz4)よくわからないのでマニュアルでノーマライズした方がいいかも。

```
usage: lefse_format_input.py [-h] [--output_table OUTPUT_TABLE] [-f {c,r}] [-c [1..n_feats]]
                             [-s [1..n_feats]] [-o float] [-u [1..n_feats]] [-m {f,s}] [-n int]
                             [-biom_c BIOM_CLASS] [-biom_s BIOM_SUBCLASS]
                             INPUT_FILE OUTPUT_FILE

LEfSe formatting modules

positional arguments:
  INPUT_FILE            the input file, feature hierarchical level can be specified with | or .
                        and those symbols must not be present for other reasons in the input file.
  OUTPUT_FILE           the output file containing the data for LEfSe

options:
  -h, --help            show this help message and exit
  --output_table OUTPUT_TABLE
                        the formatted table in txt format
  -f {c,r}              set whether the features are on rows (default) or on columns
  -c [1..n_feats]       set which feature use as class (default 1)
  -s [1..n_feats]       set which feature use as subclass (default -1 meaning no subclass)
  -o float              set the normalization value (default -1.0 meaning no normalization)
  -u [1..n_feats]       set which feature use as subject (default -1 meaning no subject)
  -m {f,s}              set the policy to adopt with missing values: f removes the features with
                        missing values, s removes samples with missing values (default f)
  -n int                set the minimum cardinality of each subclass (subclasses with low
                        cardinalities will be grouped together, if the cardinality is still low,
                        no pairwise comparison will be performed with them)
  -biom_c BIOM_CLASS    For biom input files: Set which feature use as class
  -biom_s BIOM_SUBCLASS
                        For biom input files: set which feature use as subclass
```


### lefse_run.py -h


```
usage: lefse_run.py [-h] [-o str] [-a float] [-w float] [-l float]
                    [--nlogs int] [--verbose int] [--wilc int] [-r str]
                    [--svm_norm int] [-b int] [-e int] [-c int] [-f float]
                    [-s {0,1,2}] [--min_c int] [-t str] [-y {0,1}]
                    INPUT_FILE OUTPUT_FILE

LEfSe 1.1.01

positional arguments:
  INPUT_FILE      the input file
  OUTPUT_FILE     the output file containing the data for the visualization
                  module

optional arguments:
  -h, --help      show this help message and exit
  -o str          set the file for exporting the result (only concise textual
                  form)
  -a float        set the alpha value for the Anova test (default 0.05)
  -w float        set the alpha value for the Wilcoxon test (default 0.05)
  -l float        set the threshold on the absolute value of the logarithmic
                  LDA score (default 2.0)
  --nlogs int     max log ingluence of LDA coeff
  --verbose int   verbose execution (default 0)
  --wilc int      wheter to perform the Wicoxon step (default 1)
  -r str          select LDA or SVM for effect size (default LDA)
  --svm_norm int  whether to normalize the data in [0,1] for SVM feature
                  waiting (default 1 strongly suggested)
  -b int          set the number of bootstrap iteration for LDA (default 30)
  -e int          set whether perform the wilcoxon test only among the
                  subclasses with the same name (default 0)
  -c int          set whether perform the wilcoxon test ing the Curtis's
                  approach [BETA VERSION] (default 0)
  -f float        set the subsampling fraction value for each bootstrap
                  iteration (default 0.66666)
  -s {0,1,2}      set the multiple testing correction options. 0 no correction
                  (more strict, default), 1 correction for independent
                  comparisons, 2 correction for dependent comparison
  --min_c int     minimum number of samples per subclass for performing
                  wilcoxon test (default 10)
  -t str          set the title of the analysis (default input file without
                  extension)
  -y {0,1}        (for multiclass tasks) set whether the test is performed in
                  a one-against-one ( 1 - more strict!) or in a one-against-
                  all setting ( 0 - less strict) (default 0)
```


### lefse_plot_res.py -h


```
usage: lefse_plot_res.py [-h] [--feature_font_size FEATURE_FONT_SIZE]
                         [--format {png,svg,pdf}] [--dpi DPI] [--title TITLE]
                         [--title_font_size TITLE_FONT_SIZE]
                         [--class_legend_font_size CLASS_LEGEND_FONT_SIZE]
                         [--width WIDTH] [--height HEIGHT] [--left_space LS]
                         [--right_space RS] [--orientation {h,v}]
                         [--autoscale {0,1}] [--background_color {k,w}]
                         [--subclades N_SCL]
                         [--max_feature_len MAX_FEATURE_LEN]
                         [--all_feats ALL_FEATS] [--otu_only]
                         [--report_features]
                         INPUT_FILE OUTPUT_FILE

Plot results

positional arguments:
  INPUT_FILE            tab delimited input file
  OUTPUT_FILE           the file for the output image

optional arguments:
  -h, --help            show this help message and exit
  --feature_font_size FEATURE_FONT_SIZE
                        the file for the output image
  --format {png,svg,pdf}
                        the format for the output file
  --dpi DPI
  --title TITLE
  --title_font_size TITLE_FONT_SIZE
  --class_legend_font_size CLASS_LEGEND_FONT_SIZE
  --width WIDTH
  --height HEIGHT       only for vertical histograms
  --left_space LS
  --right_space RS
  --orientation {h,v}
  --autoscale {0,1}
  --background_color {k,w}
                        set the color of the background
  --subclades N_SCL     number of label levels to be dislayed (starting from
                        the leaves, -1 means all the levels, 1 is default )
  --max_feature_len MAX_FEATURE_LEN
                        Maximum length of feature strings (def 60)
  --all_feats ALL_FEATS
  --otu_only            Plot only species resolved OTUs (as opposed to all
                        levels)
  --report_features     Report important features to STDOUT
```

### lefse_plot_cladogram.py -h

```
usage: lefse_plot_cladogram.py [-h] [--clade_sep CLADE_SEP]
                               [--max_lev MAX_LEV]
                               [--max_point_size MAX_POINT_SIZE]
                               [--min_point_size MIN_POINT_SIZE]
                               [--point_edge_width MARKEREDGEWIDTH]
                               [--siblings_connector_width SIBLINGS_CONNECTOR_WIDTH]
                               [--parents_connector_width PARENTS_CONNECTOR_WIDTH]
                               [--radial_start_lev RADIAL_START_LEV]
                               [--labeled_start_lev LABELED_START_LEV]
                               [--labeled_stop_lev LABELED_STOP_LEV]
                               [--abrv_start_lev ABRV_START_LEV]
                               [--abrv_stop_lev ABRV_STOP_LEV]
                               [--expand_void_lev EXPAND_VOID_LEV]
                               [--class_legend_vis CLASS_LEGEND_VIS]
                               [--colored_connector COLORED_CONNECTORS]
                               [--alpha ALPHA] [--title TITLE]
                               [--sub_clade SUB_CLADE]
                               [--title_font_size TITLE_FONT_SIZE]
                               [--right_space_prop R_PROP]
                               [--left_space_prop L_PROP]
                               [--label_font_size LABEL_FONT_SIZE]
                               [--background_color {k,w}]
                               [--colored_labels {0,1}]
                               [--class_legend_font_size CLASS_LEGEND_FONT_SIZE]
                               [--dpi DPI] [--format {png,svg,pdf}]
                               [--all_feats ALL_FEATS]
                               INPUT_FILE OUTPUT_FILE

Cladoplot

positional arguments:
  INPUT_FILE            tab delimited input file
  OUTPUT_FILE           the file for the output image

optional arguments:
  -h, --help            show this help message and exit
  --clade_sep CLADE_SEP
  --max_lev MAX_LEV
  --max_point_size MAX_POINT_SIZE
  --min_point_size MIN_POINT_SIZE
  --point_edge_width MARKEREDGEWIDTH
  --siblings_connector_width SIBLINGS_CONNECTOR_WIDTH
  --parents_connector_width PARENTS_CONNECTOR_WIDTH
  --radial_start_lev RADIAL_START_LEV
  --labeled_start_lev LABELED_START_LEV
  --labeled_stop_lev LABELED_STOP_LEV
  --abrv_start_lev ABRV_START_LEV
  --abrv_stop_lev ABRV_STOP_LEV
  --expand_void_lev EXPAND_VOID_LEV
  --class_legend_vis CLASS_LEGEND_VIS
  --colored_connector COLORED_CONNECTORS
  --alpha ALPHA
  --title TITLE
  --sub_clade SUB_CLADE
  --title_font_size TITLE_FONT_SIZE
  --right_space_prop R_PROP
  --left_space_prop L_PROP
  --label_font_size LABEL_FONT_SIZE
  --background_color {k,w}
                        set the color of the background
  --colored_labels {0,1}
                        draw the label with class color (1) or in black (0)
  --class_legend_font_size CLASS_LEGEND_FONT_SIZE
  --dpi DPI
  --format {png,svg,pdf}
                        the format for the output file
  --all_feats ALL_FEATS
```

```


## Example

lefse_run.pyの出力ファイルは、feature | log of the highest class average | highest class | LDA effect size | p-value
```
# 1) タブ区切りの入力ファイルから読み込み用バイナリを作成（ファイルサイズは元のタブ区切りと同じくらい）
lefse_format_input.py example.txt example.in -c 1 -s 2 -u 3 -o 1000000

# 2) LEfSeの実行
lefse_run.py -t title example.in example.res
sbatch -c 8 -J lefse_test --wrap="lefse_run.py -t example1 example.in example.res"

# 3) グラフ作成
lefse_plot_res.py       --format pdf example.res example.pdf
lefse_plot_cladogram.py --format pdf --title title --class_legend_font_size 6 example.res example_cladogram.pdf 
```# LEfSe


## Install

[LEfSe](https://huttenhower.sph.harvard.edu/lefse/)
[LEfSeの使用法](https://note.com/ytomy/n/n39b64afd3862)

[condaでLEfSeを導入するときに注意すること](https://qiita.com/FlaskSakaguchi/items/b9e70129793174d939e1)

```
conda create -n lefse python=3.7 -y
conda activate lefse
conda install -c bioconda -c conda-forge lefse -y
```


`conda create -n lefse -c bioconda lefse`だと`lefse_run.py`でエラーが出た（pythonのバージョンの問題？）

```
lefse_run.py
# Traceback (most recent call last):
#   File "/user/ngsdata/kmatsumoto/miniconda3/envs/lefse/lib/python3.10/site-packages/rpy2/rinterface_lib/openrlib.py", line 17, in <module>
#     import _rinterface_cffi_api as _rinterface_cffi  # type: ignore
# ImportError: libffi.so.7: cannot open shared object file: No such file or directory
```


## Input format

> The text tab-delimited input file consists of a list of numerical features, the class vector and optionally the subclass and subject vectors. The features can be read counts directly or abundance floating-point values more generally, and the first field is the name of the feature. Class, subclass and subject vectors have a name (the first field) and a list of non-numerical strings.

> Although both column and row feature organization is accepted, given the high-dimensional nature of metagenomic data, the listing of the features in rows is preferred.

- タブ区切りの数値表。
- １列目はラベル（文字列）。
- １行目がグループ名、２行目？がサブグループ名（任意）、３行目がID（これらも１列目に名称をつける）。４行目以降が特徴量。
- 数値はカウントでも頻度でも。
- 行が分類群や遺伝子で、列がサンプル。
- 分類群名は"|"で区切ると上の分類群も検定してくれる。（上の分類群の行がなくても勝手に集計してくれるもよう。行を入れる場合は合計頻度を入れるもよう（サンプルファイル参照）。）。"|"で区切るとcladogramも出せる。これらをしてほしくない時は|で区切らなければOK。
- 比較は１行目のグループ間で行われる。サブグループは、サブグループごとでも（それぞれのサブグループ内でも）グループ間で違うか、というように使われるみたい。[LefSe Class and Subclass](https://groups.google.com/g/lefse-users/c/QPT7S-o-RhM)


例）[sample input file: hmp_small_aerobiosis.txt](https://github.com/biobakery/biobakery/raw/master/test_suite/biobakery_tests/data/lefse/input/hmp_small_aerobiosis.txt)
```
oxygen_availability	High_O2	Mid_O2	Low_O2	Mid_O2	Low_O2
body_site	ear	oral	gut	oral	gut
subject_id	158721788	158721788	159146620	159005010	159166850
Archaea|Euryarchaeota|Methanobacteria|Methanobacteriales|Methanobacteriaceae|Methanobrevibacter	2.96541e-06	5.08937e-06	4.93921e-06	7.86782e-06	6.22604e-06
Bacteria	0.999994	0.99999	0.99999	0.999984	0.999988
Bacteria|Acidobacteria	5.0412e-05	8.65194e-05	8.39666e-05	0.000133753	0.000105843
Bacteria|Acidobacteria|Acidobacteria_Gp10|Gp10	2.96541e-06	5.08937e-06	4.93921e-06	7.86782e-06	6.22604e-06
Bacteria|Acidobacteria|Acidobacteria_Gp11|Gp11	2.96541e-06	5.08937e-06	4.93921e-06	7.86782e-06	6.22604e-06
Bacteria|Acidobacteria|Acidobacteria_Gp16|Gp16	2.96541e-06	5.08937e-06	4.93921e-06	7.86782e-06	6.22604e-06
Bacteria|Acidobacteria|Acidobacteria_Gp17|Gp17	2.96541e-06	5.08937e-06	4.93921e-06	7.86782e-06	6.22604e-06
```


## ヘルプメッセージ

### lefse_format_input.py -h

- normalizationは、サンプル（列）ごとに列の合計値で割った後に`-o`の値を掛ける？[Questions about pre-sample normalization and low aboundant taxa for LEfSe](https://groups.google.com/g/lefse-users/c/5IFgZRjHGz4)よくわからないのでマニュアルでノーマライズした方がいいかも。

```
usage: lefse_format_input.py [-h] [--output_table OUTPUT_TABLE] [-f {c,r}] [-c [1..n_feats]]
                             [-s [1..n_feats]] [-o float] [-u [1..n_feats]] [-m {f,s}] [-n int]
                             [-biom_c BIOM_CLASS] [-biom_s BIOM_SUBCLASS]
                             INPUT_FILE OUTPUT_FILE

LEfSe formatting modules

positional arguments:
  INPUT_FILE            the input file, feature hierarchical level can be specified with | or .
                        and those symbols must not be present for other reasons in the input file.
  OUTPUT_FILE           the output file containing the data for LEfSe

options:
  -h, --help            show this help message and exit
  --output_table OUTPUT_TABLE
                        the formatted table in txt format
  -f {c,r}              set whether the features are on rows (default) or on columns
  -c [1..n_feats]       set which feature use as class (default 1)
  -s [1..n_feats]       set which feature use as subclass (default -1 meaning no subclass)
  -o float              set the normalization value (default -1.0 meaning no normalization)
  -u [1..n_feats]       set which feature use as subject (default -1 meaning no subject)
  -m {f,s}              set the policy to adopt with missing values: f removes the features with
                        missing values, s removes samples with missing values (default f)
  -n int                set the minimum cardinality of each subclass (subclasses with low
                        cardinalities will be grouped together, if the cardinality is still low,
                        no pairwise comparison will be performed with them)
  -biom_c BIOM_CLASS    For biom input files: Set which feature use as class
  -biom_s BIOM_SUBCLASS
                        For biom input files: set which feature use as subclass
```


### lefse_run.py -h


```
usage: lefse_run.py [-h] [-o str] [-a float] [-w float] [-l float]
                    [--nlogs int] [--verbose int] [--wilc int] [-r str]
                    [--svm_norm int] [-b int] [-e int] [-c int] [-f float]
                    [-s {0,1,2}] [--min_c int] [-t str] [-y {0,1}]
                    INPUT_FILE OUTPUT_FILE

LEfSe 1.1.01

positional arguments:
  INPUT_FILE      the input file
  OUTPUT_FILE     the output file containing the data for the visualization
                  module

optional arguments:
  -h, --help      show this help message and exit
  -o str          set the file for exporting the result (only concise textual
                  form)
  -a float        set the alpha value for the Anova test (default 0.05)
  -w float        set the alpha value for the Wilcoxon test (default 0.05)
  -l float        set the threshold on the absolute value of the logarithmic
                  LDA score (default 2.0)
  --nlogs int     max log ingluence of LDA coeff
  --verbose int   verbose execution (default 0)
  --wilc int      wheter to perform the Wicoxon step (default 1)
  -r str          select LDA or SVM for effect size (default LDA)
  --svm_norm int  whether to normalize the data in [0,1] for SVM feature
                  waiting (default 1 strongly suggested)
  -b int          set the number of bootstrap iteration for LDA (default 30)
  -e int          set whether perform the wilcoxon test only among the
                  subclasses with the same name (default 0)
  -c int          set whether perform the wilcoxon test ing the Curtis's
                  approach [BETA VERSION] (default 0)
  -f float        set the subsampling fraction value for each bootstrap
                  iteration (default 0.66666)
  -s {0,1,2}      set the multiple testing correction options. 0 no correction
                  (more strict, default), 1 correction for independent
                  comparisons, 2 correction for dependent comparison
  --min_c int     minimum number of samples per subclass for performing
                  wilcoxon test (default 10)
  -t str          set the title of the analysis (default input file without
                  extension)
  -y {0,1}        (for multiclass tasks) set whether the test is performed in
                  a one-against-one ( 1 - more strict!) or in a one-against-
                  all setting ( 0 - less strict) (default 0)
```


### lefse_plot_res.py -h


```
usage: lefse_plot_res.py [-h] [--feature_font_size FEATURE_FONT_SIZE]
                         [--format {png,svg,pdf}] [--dpi DPI] [--title TITLE]
                         [--title_font_size TITLE_FONT_SIZE]
                         [--class_legend_font_size CLASS_LEGEND_FONT_SIZE]
                         [--width WIDTH] [--height HEIGHT] [--left_space LS]
                         [--right_space RS] [--orientation {h,v}]
                         [--autoscale {0,1}] [--background_color {k,w}]
                         [--subclades N_SCL]
                         [--max_feature_len MAX_FEATURE_LEN]
                         [--all_feats ALL_FEATS] [--otu_only]
                         [--report_features]
                         INPUT_FILE OUTPUT_FILE

Plot results

positional arguments:
  INPUT_FILE            tab delimited input file
  OUTPUT_FILE           the file for the output image

optional arguments:
  -h, --help            show this help message and exit
  --feature_font_size FEATURE_FONT_SIZE
                        the file for the output image
  --format {png,svg,pdf}
                        the format for the output file
  --dpi DPI
  --title TITLE
  --title_font_size TITLE_FONT_SIZE
  --class_legend_font_size CLASS_LEGEND_FONT_SIZE
  --width WIDTH
  --height HEIGHT       only for vertical histograms
  --left_space LS
  --right_space RS
  --orientation {h,v}
  --autoscale {0,1}
  --background_color {k,w}
                        set the color of the background
  --subclades N_SCL     number of label levels to be dislayed (starting from
                        the leaves, -1 means all the levels, 1 is default )
  --max_feature_len MAX_FEATURE_LEN
                        Maximum length of feature strings (def 60)
  --all_feats ALL_FEATS
  --otu_only            Plot only species resolved OTUs (as opposed to all
                        levels)
  --report_features     Report important features to STDOUT
```

### lefse_plot_cladogram.py -h

```
usage: lefse_plot_cladogram.py [-h] [--clade_sep CLADE_SEP]
                               [--max_lev MAX_LEV]
                               [--max_point_size MAX_POINT_SIZE]
                               [--min_point_size MIN_POINT_SIZE]
                               [--point_edge_width MARKEREDGEWIDTH]
                               [--siblings_connector_width SIBLINGS_CONNECTOR_WIDTH]
                               [--parents_connector_width PARENTS_CONNECTOR_WIDTH]
                               [--radial_start_lev RADIAL_START_LEV]
                               [--labeled_start_lev LABELED_START_LEV]
                               [--labeled_stop_lev LABELED_STOP_LEV]
                               [--abrv_start_lev ABRV_START_LEV]
                               [--abrv_stop_lev ABRV_STOP_LEV]
                               [--expand_void_lev EXPAND_VOID_LEV]
                               [--class_legend_vis CLASS_LEGEND_VIS]
                               [--colored_connector COLORED_CONNECTORS]
                               [--alpha ALPHA] [--title TITLE]
                               [--sub_clade SUB_CLADE]
                               [--title_font_size TITLE_FONT_SIZE]
                               [--right_space_prop R_PROP]
                               [--left_space_prop L_PROP]
                               [--label_font_size LABEL_FONT_SIZE]
                               [--background_color {k,w}]
                               [--colored_labels {0,1}]
                               [--class_legend_font_size CLASS_LEGEND_FONT_SIZE]
                               [--dpi DPI] [--format {png,svg,pdf}]
                               [--all_feats ALL_FEATS]
                               INPUT_FILE OUTPUT_FILE

Cladoplot

positional arguments:
  INPUT_FILE            tab delimited input file
  OUTPUT_FILE           the file for the output image

optional arguments:
  -h, --help            show this help message and exit
  --clade_sep CLADE_SEP
  --max_lev MAX_LEV
  --max_point_size MAX_POINT_SIZE
  --min_point_size MIN_POINT_SIZE
  --point_edge_width MARKEREDGEWIDTH
  --siblings_connector_width SIBLINGS_CONNECTOR_WIDTH
  --parents_connector_width PARENTS_CONNECTOR_WIDTH
  --radial_start_lev RADIAL_START_LEV
  --labeled_start_lev LABELED_START_LEV
  --labeled_stop_lev LABELED_STOP_LEV
  --abrv_start_lev ABRV_START_LEV
  --abrv_stop_lev ABRV_STOP_LEV
  --expand_void_lev EXPAND_VOID_LEV
  --class_legend_vis CLASS_LEGEND_VIS
  --colored_connector COLORED_CONNECTORS
  --alpha ALPHA
  --title TITLE
  --sub_clade SUB_CLADE
  --title_font_size TITLE_FONT_SIZE
  --right_space_prop R_PROP
  --left_space_prop L_PROP
  --label_font_size LABEL_FONT_SIZE
  --background_color {k,w}
                        set the color of the background
  --colored_labels {0,1}
                        draw the label with class color (1) or in black (0)
  --class_legend_font_size CLASS_LEGEND_FONT_SIZE
  --dpi DPI
  --format {png,svg,pdf}
                        the format for the output file
  --all_feats ALL_FEATS
```

```


## Example

lefse_run.pyの出力ファイルは、feature | log of the highest class average | highest class | LDA effect size | p-value
```
# 1) タブ区切りの入力ファイルから読み込み用バイナリを作成（ファイルサイズは元のタブ区切りと同じくらい）
lefse_format_input.py example.txt example.in -c 1 -s 2 -u 3 -o 1000000

# 2) LEfSeの実行
lefse_run.py -t title example.in example.res
sbatch -c 8 -J lefse_test --wrap="lefse_run.py -t example1 example.in example.res"

# 3) グラフ作成
lefse_plot_res.py       --format pdf example.res example.pdf
lefse_plot_cladogram.py --format pdf --title title --class_legend_font_size 6 example.res example_cladogram.pdf 
```