# CheckM

[CheckM](https://ecogenomics.github.io/CheckM/)
[Ecogenomics/CheckM/wiki](https://github.com/Ecogenomics/CheckM/wiki)

```
conda create -n checkm python=3.9 -y
conda activate checkm
conda install numpy matplotlib pysam -y
conda install hmmer prodigal pplacer -y
conda install checkm-genome -y
```

condaでインストールした場合、DBも用意してくれる。場所は~/miniconda3/envs/checkm/checkm_data。


```
>checkm lineage_wf -h
usage: checkm lineage_wf [-h] [-r] [--ali] [--nt] [-g] [-u UNIQUE] [-m MULTI] [--force_domain]
                         [--no_refinement] [--individual_markers] [--skip_adj_correction]
                         [--skip_pseudogene_correction] [--aai_strain AAI_STRAIN]
                         [-a ALIGNMENT_FILE] [--ignore_thresholds] [-e E_VALUE] [-l LENGTH]
                         [-f FILE] [--tab_table] [-x EXTENSION] [-t THREADS]
                         [--pplacer_threads PPLACER_THREADS] [-q] [--tmpdir TMPDIR]
                         bin_input output_dir

Runs tree, lineage_set, analyze, qa

positional arguments:
  bin_input             directory containing bins (fasta format) or path to file describing genomes/genes - tab separated in 2 or 3 columns [genome ID, genome fna, genome translation file (pep)]
  output_dir            directory to write output files

optional arguments:
  -h, --help            show this help message and exit
  -r, --reduced_tree    use reduced tree (requires <16GB of memory) for determining lineage of each bin
  --ali                 generate HMMER alignment file for each bin
  --nt                  generate nucleotide gene sequences for each bin
  -g, --genes           bins contain genes as amino acids instead of nucleotide contigs
  -u, --unique UNIQUE   minimum number of unique phylogenetic markers required to use lineage-specific marker set (default: 10)
  -m, --multi MULTI     maximum number of multi-copy phylogenetic markers before defaulting to domain-level marker set (default: 10)
  --force_domain        use domain-level sets for all bins
  --no_refinement       do not perform lineage-specific marker set refinement
  --individual_markers  treat marker as independent (i.e., ignore co-located set structure)
  --skip_adj_correction
                        do not exclude adjacent marker genes when estimating contamination
  --skip_pseudogene_correction
                        skip identification and filtering of pseudogenes
  --aai_strain AAI_STRAIN
                        AAI threshold used to identify strain heterogeneity (default: 0.9)
  -a, --alignment_file ALIGNMENT_FILE
                        produce file showing alignment of multi-copy genes and their AAI identity
  --ignore_thresholds   ignore model-specific score thresholds
  -e, --e_value E_VALUE
                        e-value cut off (default: 1e-10)
  -l, --length LENGTH   percent overlap between target and query (default: 0.7)
  -f, --file FILE       print results to file (default: stdout)
  --tab_table           print tab-separated values table
  -x, --extension EXTENSION
                        extension of bins (other files in directory are ignored) (default: fna)
  -t, --threads THREADS
                        number of threads (default: 1)
  --pplacer_threads PPLACER_THREADS
                        number of threads used by pplacer (memory usage increases linearly with additional threads) (default: 1)
  -q, --quiet           suppress console output
  --tmpdir TMPDIR       specify an alternative directory for temporary files

Example: checkm lineage_wf ./bins ./output
```

```
sbatch -c 22 -p c6420 --wrap="checkm lineage_wf \
 -t 22 \
 --pplacer_threads 22 \
 --extension fasta /imetgpfs/projects/cw/KOEI/220413-No1368-n294/02_abricate/test/ test_out | tee log.txt"

```