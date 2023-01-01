
```
get_blast_tophit <- function(blast_file){
	command <- paste("cat", blast_file, "|", "awk -F '\t' '!x[$1]++'")
	blast_tophit <- read.table(text = system(command, intern=TRUE), sep = "\t", quote="") %>%
		as_tibble() %>%
		rename(paste0("V", 1:12) %>% set_names(c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen","qstart", "qend", "sstart", "send", "evalue", "bitscore")))
	return(blast_tophit)
}
```