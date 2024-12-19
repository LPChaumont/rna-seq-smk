log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library("tximport")
library("readr")
library("jsonlite")


write_tsv <- function(data, file) {
  write.table(
    data, file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

add_gene_col <- function(df) {
  return(data.frame(gene_id = rownames(df), df))
}

add_tx_col <- function(df) {
  return(data.frame(tx_id = gsub("\\.[0-9]*$", "", rownames(df)), df))
}


quant_files <- snakemake@input[["quants"]]
sample_names <- sapply(quant_files, \(x) basename(dirname(x)))
names(quant_files) <- sample_names

tx2gene <- read.table(
  snakemake@input[["tx2gene"]],
  header= TRUE,
  sep="\t",
  check.names = FALSE
)

txi_tx <- tximport(
  quant_files,
  type = "salmon",
  tx2gene = tx2gene,
  ignoreTxVersion = TRUE,
  txOut = TRUE
)
saveRDS(txi_tx, file=snakemake@output[["tx_rds"]])

txi_gene <- tximport(
  quant_files,
  type = "salmon",
  tx2gene = tx2gene,
  ignoreTxVersion = TRUE,
  txOut = FALSE
)
saveRDS(txi_gene, file=snakemake@output[["gene_rds"]])

write_tsv(add_tx_col(txi_tx$counts), snakemake@output[["tx_counts"]])
write_tsv(add_gene_col(txi_gene$counts), snakemake@output[["gene_counts"]])
write_tsv(add_tx_col(txi_tx$abundance), snakemake@output[["tx_tpm"]])
write_tsv(add_gene_col(txi_gene$abundance), snakemake@output[["gene_tpm"]])
