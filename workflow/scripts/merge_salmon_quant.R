log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library("tximport")
library("readr")
library("jsonlite")


write_tsv <- function(data, file) {
  write.table(
    data.frame(gene_id = row.names(data), data),
    file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}


quant_files <- list.files(
  snakemake@params[["quantdir"]],
  pattern = "quant.sf$",
  recursive = TRUE,
  full.names = TRUE
)
sample_names <- basename(dirname(quant_files))
names(quant_files) <- sample_names

tx2gene <- read.table(
  snakemake@input[["tx2gene"]],
  header= TRUE,
  sep="\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

txi_tx <- tximport(
  quant_files,
  type = "salmon",
  tx2gene = tx2gene,
  ignoreTxVersion = TRUE
  txOut = TRUE
)
saveRDS(txi_tx, file=snakemake@output[["tx_rds"]])

txi_gene <- tximport(
  quant_files,
  type = "salmon",
  tx2gene = tx2gene,
  txOut = FALSE
)
saveRDS(txi_gene, file=snakemake@output[["gene_rds"]])

write_tsv(txi_tx$counts, snakemake@output[["tx_count"]])
write_tsv(txi_gene$counts, snakemake@output[["gene_count"]])
write_tsv(txi_tx$abundance, snakemake@output[["tx_tpm"]])
write_tsv(txi_gene$abundance, snakemake@output[["gene_tpm"]])
