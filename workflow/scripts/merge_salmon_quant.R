log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(tximport)
library(rtracklayer)

# Get salmon quant files
quants_dir <- snakemake@params[["quantdir"]]
quant_files <- list.files(quants_dir, pattern = "quant.sf$", recursive = TRUE, full.names = TRUE)
sample_names <- basename(dirname(quant_files))
names(quant_files) <- sample_names

# Transcript id mapping to gene id
gtf <- rtracklayer::import(snakemake@input[["gtf"]])
tx2gene <- mcols(gtf[gtf$type == "transcript"][, c("transcript_id", "gene_id")])

# Make tx2gene with GenomicFeatures
# library(GenomicFeatures)
# txdb <- makeTxDbFromGFF(gtf_path, format = "gtf")
# k <- keys(txdb, keytype = "GENEID")
# tx2gene <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")

# tximport
transcript <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
gene <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, txOut = TRUE)


write_tsv <- function(data, filename) {
  write.table(data.frame(gene_id = row.names(data), data),
    filename,
    sep = "\t", row.names = FALSE, quote = FALSE
  )
}

write_tsv(transcript$counts, snakemake@output[["salmon_transcript_counts.tsv"]])
write_tsv(transcript$abundance, snakemake@output[["salmon_transcript_tpm.tsv"]])
write_tsv(gene$counts, snakemake@output[["salmon_gene_counts.tsv"]])
write_tsv(gene$abundance, snakemake@output[["salmon_gene_tpm.tsv"]])
