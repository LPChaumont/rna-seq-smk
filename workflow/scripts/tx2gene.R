log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(GenomicFeatures)

txdb <- makeTxDbFromGFF(snakemake@input[[1]], format = "gtf")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

write.table(tx2gene, file = snakemake@output[[1]], quote = FALSE, row.names = FALSE, sep = "\t")
