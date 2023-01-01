# argparser

依存が少ないのでRの引数受け取りはこちらがおすすめ。

```
	parser <- arg_parser(
		"Conduct edgeR.",
		hide.opts = TRUE
	)
	parser <- parser %>%
		add_argument("--infile",    "REQUIRED. File name of count matrix. 1st column contains gene names, and subsequent columns are count data.") %>%
		add_argument("--group",     "Group names of the samples (e.g. --group treat1 treat1 treat2 treat2 ctl ctl)", nargs=Inf) %>%
		add_argument("--suffix",    "When --group is not specified, group names are automatically obtained by removing this suffix from the column names.", default="_[0-9]+$") %>%
		add_argument("--compare",   "Hyphen separated group names you want to compare. (e.g. --compare treat1-ctl treat2-ctl). If not specified, each combination will be compared.", nargs=Inf) %>%
		add_argument("--min_count", "Minimum count required for at least some samples. Low count genes will not be used for testing. (see edger::filterByExpr).", default=10) %>%
		add_argument("--test",      "[QLF, LRT]; quasi-likelihood F-test or likelihood ratio test.", default="QLF") %>%
		add_argument("--drop",      "Rows with this value in the 1st column will be dropped at the begining. (e.g. Unannotated)", nargs=Inf)
	args <- parse_args(parser)

	cat("\n==== edger.R ====\n")

	stopifnot("Args --infile is required." = !is.na(args$infile))
	if(!file.exists(args$infile)){stop(paste0("'", args$infile, "' does not exist."))}

	# Read & modify count data file
	# read.table instead of read_tsv in preparation for data file with index (read_tsv raise warning because of missing column name))
	count_matrix <- read.table(args$infile, sep = "\t", header = TRUE, row.names = 1)
	stopifnot("The second and subsequent columns should be numeric." = (all(map_lgl(count_matrix, ~ is.numeric(.)))))
	if(all(!is.na(args$drop))){
		count_matrix <- count_matrix[!(rownames(count_matrix) %in% args$drop),]
	}

	res <- edger(count = count_matrix, group = args$group, suffix = args$suffix, compare = args$compare, min_count = args$min_count, test = args$test)

```

コマンドから呼び出した時だけ実行する場合は、`if(sys.nframe() == 0){}`の中に入れる。
引数のデフォルトは、`defailt=`で指定しなかった場合NAが入る。
