# ncbi-genome-download

[kblin/ncbi-genome-download](https://github.com/kblin/ncbi-genome-download)
[NCBI FTPサーバからゲノム配列をダウンロードする ncbi-genome-download](https://kazumaxneo.hatenablog.com/entry/2018/12/08/073000)

```
conda install -c bioconda ncbi-genome-download
```


ヘルプ

```
usage: ncbi-genome-download [-h] [-s {refseq,genbank}] [-F FILE_FORMATS] [-l ASSEMBLY_LEVELS]
                            [-g GENERA] [--genus GENERA] [--fuzzy-genus] [-S STRAINS]
                            [-T SPECIES_TAXIDS] [-t TAXIDS] [-A ASSEMBLY_ACCESSIONS]
                            [-R REFSEQ_CATEGORIES] [--refseq-category REFSEQ_CATEGORIES]
                            [-o OUTPUT] [--flat-output] [-H] [-P] [-u URI] [-p N] [-r N]
                            [-m METADATA_TABLE] [-n] [-N] [-v] [-d] [-V] [-M TYPE_MATERIALS]
                            groups

positional arguments:
  groups                The NCBI taxonomic groups to download (default: all). A comma-separated
                        list of taxonomic groups is also possible. For example:
                        "bacteria,viral"Choose from: ['all', 'archaea', 'bacteria', 'fungi',
                        'invertebrate', 'metagenomes', 'plant', 'protozoa',
                        'vertebrate_mammalian', 'vertebrate_other', 'viral']

optional arguments:
  -h, --help            show this help message and exit
  -s {refseq,genbank}, --section {refseq,genbank}
                        NCBI section to download (default: refseq)
  -F FILE_FORMATS, --formats FILE_FORMATS
                        Which formats to download (default: genbank).A comma-separated list of
                        formats is also possible. For example: "fasta,assembly-report". Choose
                        from: ['genbank', 'fasta', 'rm', 'features', 'gff', 'protein-fasta',
                        'genpept', 'wgs', 'cds-fasta', 'rna-fna', 'rna-fasta', 'assembly-report',
                        'assembly-stats', 'all']
  -l ASSEMBLY_LEVELS, --assembly-levels ASSEMBLY_LEVELS
                        Assembly levels of genomes to download (default: all). A comma-separated
                        list of assembly levels is also possible. For example:
                        "complete,chromosome". Choose from: ['all', 'complete', 'chromosome',
                        'scaffold', 'contig']
  -g GENERA, --genera GENERA
                        Only download sequences of the provided genera. A comma-seperated list of
                        genera is also possible. For example: "Streptomyces coelicolor,Escherichia
                        coli". (default: [])
  --genus GENERA        Deprecated alias of --genera
  --fuzzy-genus         Use a fuzzy search on the organism name instead of an exact match.
  -S STRAINS, --strains STRAINS
                        Only download sequences of the given strain(s). A comma-separated list of
                        strain names is possible, as well as a path to a filename containing one
                        name per line.
  -T SPECIES_TAXIDS, --species-taxids SPECIES_TAXIDS
                        Only download sequences of the provided species NCBI taxonomy IDs. A
                        comma-separated list of species taxids is also possible. For example:
                        "52342,12325". (default: [])
  -t TAXIDS, --taxids TAXIDS
                        Only download sequences of the provided NCBI taxonomy IDs. A comma-
                        separated list of taxids is also possible. For example: "9606,9685".
                        (default: [])
  -A ASSEMBLY_ACCESSIONS, --assembly-accessions ASSEMBLY_ACCESSIONS
                        Only download sequences matching the provided NCBI assembly accession(s).
                        A comma-separated list of accessions is possible, as well as a path to a
                        filename containing one accession per line.
  -R REFSEQ_CATEGORIES, --refseq-categories REFSEQ_CATEGORIES
                        Only download sequences of the provided refseq categories (default: all)
  --refseq-category REFSEQ_CATEGORIES
                        Deprecated alias for --refseq-categories
  -o OUTPUT, --output-folder OUTPUT
                        Create output hierarchy in specified folder (default:
                        /auto/user/ngsdata/kmatsumoto/jutaku/220420_NKTA1187)
  --flat-output         Dump all files right into the output folder without creating any
                        subfolders.
  -H, --human-readable  Create links in human-readable hierarchy (might fail on Windows)
  -P, --progress-bar    Create a progress bar for indicating the download progress
  -u URI, --uri URI     NCBI base URI to use (default: https://ftp.ncbi.nih.gov/genomes)
  -p N, --parallel N    Run N downloads in parallel (default: 1)
  -r N, --retries N     Retry download N times when connection to NCBI fails (default: 0)
  -m METADATA_TABLE, --metadata-table METADATA_TABLE
                        Save tab-delimited file with genome metadata
  -n, --dry-run         Only check which files to download, don't download genome files.
  -N, --no-cache        Don't cache the assembly summary file in
                        /user/ngsdata/kmatsumoto/.cache/ncbi-genome-download.
  -v, --verbose         increase output verbosity
  -d, --debug           print debugging information
  -V, --version         print version information
  -M TYPE_MATERIALS, --type-materials TYPE_MATERIALS
                        Specifies the relation to type material for the assembly (default: any).
                        "any" will include assemblies with no relation to type material value
                        defined, "all" will download only assemblies with a defined value. A
                        comma-separated list of relatons. For example: "reference,synonym". Choose
                        from: ['any', 'all', 'type', 'reference', 'synonym', 'proxytype',
                        'neotype'] .
```