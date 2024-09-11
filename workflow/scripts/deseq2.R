log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ashr)

args <- base::commandArgs(trailingOnly = TRUE)
outdir_path <- args[1]
counts_path <- args[2]
samples_path <- args[3]
contrasts_path <- args[4]
full_model <- args[5]
reduced_model <- args[6]
tpm_path <- args[7]
control_genes_path <- args[8]
padj_thresh <- as.numeric(args[9])
lfc_thresh <- as.numeric(args[10])
min_gene_expr <- as.integer(args[11])
min_samps_gene_expr <- as.integer(args[12])


read_sample_file <- function(file_path, full_model) {
  samples_df <- read.table(file_path,
    header = TRUE, sep = "\t", quote = "",
    row.names = "sample", check.names = FALSE
  )

  samples_df[] <- data.frame(lapply(samples_df, as.factor)) # the "[]" keeps the dataframe structure

  return(samples_df)
}


read_count_file <- function(file_path) {
  counts <- read.table(file_path,
    header = TRUE, sep = "\t", quote = "",
    row.names = 1, check.names = FALSE
  )

  counts <- round(counts[, sapply(counts, is.numeric)])

  return(counts)
}


read_tpm_file <- function(file_path) {
  tpm <- read.table(file_path,
    header = TRUE, sep = "\t", quote = "",
    row.names = 1, check.names = FALSE
  )

  return(tpm)
}


read_contrast_file <- function(file_path) {
  comparisons <- read.table(file_path,
    header = TRUE, sep = "\t", quote = "",
    row.names = NULL, check.names = FALSE
  )

  return(comparisons)
}


write_tsv <- function(dataframe, outdir, filename) {
  # Write a TSV file with the first column being the row names (gene_id)
  df <- data.frame(gene_id = row.names(dataframe), dataframe)
  file <- file.path(outdir, filename)
  write.table(df,
    file = file,
    row.names = FALSE,
    sep = "\t",
    quote = FALSE
  )
}


pre_filter_dds <- function(dds, tpm_path, min_gene_expr, min_samps_gene_expr, samples_df) {
  if (tpm_path != "") {
    message(
      "\nPre-filtering gene expression with TPM:\n",
      "Minimal gene expression: ", min_gene_expr, "\n",
      "Minimal number of samples where genes should be expressed: ", min_samps_gene_expr, "\n"
    )
    tpm_df <- read_tpm_file(tpm_path)

    tpm_df <- tpm_df[, sapply(tpm_df, is.numeric)]
    keep <- rowSums(tpm_df >= min_gene_expr) >= min_samps_gene_expr
    keep <- row.names(tpm_df[keep, ])
    dds <- dds[row.names(dds) %in% keep, ]
  } else {
    message(
      "\nPre-filtering gene expression with counts:\n",
      "Minimal gene expression: ", min_gene_expr, "\n",
      "Minimal number of samples where genes should be expressed: ", min_samps_gene_expr, "\n"
    )
    keep <- rowSums(counts(dds) >= min_gene_expr) >= min_samps_gene_expr
    dds <- dds[keep, ]
  }
  return(dds)
}


estimate_size_factors <- function(dds, control_genes_path, outdir_path) {
  # With control genes
  if (control_genes_path != "") {
    control_genes <- readLines(control_genes_path)
    is_control_genes <- rownames(dds) %in% control_genes
    dds <- estimateSizeFactors(dds, controlGenes = is_control_genes)

    size_factors <- sizeFactors(dds)
    size_factors_df <- data.frame(
      sample = names(size_factors),
      size_factor = size_factors
    )
    write.table(size_factors_df,
      file = file.path(outdir_path, "size_factors_with_control_genes.tsv"),
      sep = "\t", row.names = FALSE, quote = FALSE
    )
  } else {
    # Without control genes
    dds <- estimateSizeFactors(dds)
    size_factors <- sizeFactors(dds)

    size_factors_df <- data.frame(
      sample = names(size_factors),
      size_factor = size_factors
    )
    write.table(size_factors_df,
      file = file.path(outdir_path, "size_factors.tsv"),
      sep = "\t", row.names = FALSE, quote = FALSE
    )
  }
}


filter_diffexpr <- function(df, padj_thresh = 0.05, lfc_thresh = 1) {
  df <- na.omit(df)
  if (min(df$padj) == 0) {
    df$padj[df$padj == 0] <- 1e-300
  }
  df <- df[df$padj < padj_thresh & abs(df$log2FoldChange) > lfc_thresh, ]

  return(df)
}


theme_blank <- function() {
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 16),
    plot.subtitle = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )
}

plot_pca <- function(vsd, ntop) {
  pca_data <- DESeq2::plotPCA(vsd, intgroup = "condition", ntop = ntop, returnData = TRUE)
  percent_var <- round(100 * attr(pca_data, "percentVar"))
  len_vsd <- length(vsd)
  plot_subtitle <- ifelse(ntop == len_vsd, paste("All", len_vsd, "genes"), paste("Top", ntop, "genes"))

  ggplot(pca_data, aes(PC1, PC2, color = condition)) +
    geom_point(size = 3) +
    labs(
      title = "First PCs on vst-transformed data",
      subtitle = plot_subtitle,
      x = paste0("PC1 (", percent_var[1], "%)"),
      y = paste0("PC2 (", percent_var[2], "%)"),
      color = "Condition"
    ) +
    coord_fixed() +
    theme_blank()
}

plot_pvalues <- function(diffexpr_df, subtitle = NULL) {
  ggplot(diffexpr_df, aes(x = pvalue)) +
    geom_histogram(fill = "seagreen4") +
    labs(
      title = paste("p-value Distribution"),
      subtitle = subtitle,
      x = "p-values",
      y = "Frequency"
    ) +
    theme_blank()
}

plot_heatmap_samples <- function(vsd) {
  sampleDists <- dist(t(SummarizedExperiment::assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)
  colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(255)
  pheatmap::pheatmap(sampleDistMatrix,
    clustering_distance_rows = sampleDists,
    clustering_distance_cols = sampleDists,
    col = colors,
    clustering_method = "ward.D2",
    main = paste("Sample-to-sample Euclidean distances on vst-transformed data")
  )
}


# Prepare DESeqDataSet object --------------------------------------------------
if (file.exists(outdir_path) == FALSE) {
  dir.create(outdir_path, recursive = TRUE)
}

samples_df <- read_sample_file(samples_path, full_model)

counts <- read_count_file(counts_path)

contrast_df <- read_contrast_file(contrasts_path)

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts,
  colData = samples_df,
  design = as.formula(full_model)
)

saveRDS(dds, file.path(outdir_path, "raw_dds.rds"))

dds <- pre_filter_dds(dds,
  tpm_path = tpm_path,
  min_gene_expr = min_gene_expr,
  min_samps_gene_expr = min_samps_gene_expr,
  samples_df = samples_df
)

if (control_genes_path != "") {
  message("Estimating size factors\n")
} else {
  message("Estimating size factors with control genes\n")
}

estimate_size_factors(dds, control_genes_path, outdir_path)


# Differential expression analysis with the Wald test --------------------------
message("Differential expression analysis with the Wald test")
dds_wald <- DESeq2::DESeq(dds, test = "Wald")

# Get contrasts results from Wald test
for (i in 1:nrow(contrast_df)) {
  comparison <- contrast_df$contrast[i]
  reference <- contrast_df$reference_level[i]
  test <- contrast_df$tested_level[i]
  res_wald <- DESeq2::results(dds_wald, contrast = c("condition", reference, test))

  write_tsv(
    res_wald,
    outdir_path,
    paste0("differential_genes_", comparison, ".tsv")
  )

  write_tsv(
    filter_diffexpr(res_wald, padj_thresh, lfc_thresh),
    outdir_path,
    paste0("differential_genes_filtered_", comparison, ".tsv")
  )

  message("Plotting p-value histogram")
  plot_pvalues(res_wald, paste("Wald Test:", comparison))
  pvalue_file_path <- paste0(outdir_path, "/pvalue_histogram_", comparison, ".png")
  ggsave(pvalue_file_path, dpi = 600)

  # ----------------------------------------------------------------------------
  # MA plot
  # Same results if it is Wald or LRT
  message("Plotting MA plot")
  resAshT <- DESeq2::lfcShrink(dds_wald,
    contrast = c("condition", reference, test),
    type = "ashr",
    lfcThreshold = 1
  )

  png(
    filename = paste0(outdir_path, "/MA_plot_", comparison, ".png"),
    width = 7, height = 7, units = "in", res = 600
  )
  plotMA(resAshT, main = comparison, alpha = 0.05)
  abline(h = c(-1, 1), col = "dodgerblue", lwd = 2)
  dev.off()
}

# Differential expression analysis with the LRT test ---------------------------
if (reduced_model != "") {
  message("Differential expression analysis with the Likelihood Ratio Test")
  dds_lrt <- DESeq2::DESeq(dds, test = "LRT", reduced = as.formula(reduced_model))
  res_lrt <- DESeq2::results(dds_lrt)

  write_tsv(
    res_lrt,
    outdir_path,
    "differential_genes_LRT.tsv"
  )

  write_tsv(
    filter_diffexpr(res_lrt, padj_thresh, lfc_thresh),
    outdir_path,
    "differential_genes_LRT_filtered.tsv"
  )

  message("Plotting p-value histogram")
  plot_pvalues(res_lrt, "LRT Test")
  pvalue_file_path <- file.path(outdir_path, "pvalue_histogram_LRT.png")
  ggsave(pvalue_file_path, dpi = 600)
}

# ------------------------------------------------------------------------------
# Dispersion plot
# Same results if it is Wald or LRT
message("Plotting per-gene dispersion estimates")
png(
  filename = file.path(outdir_path, "dispersion_plot.png"),
  width = 7, height = 7, units = "in", res = 600
)
DESeq2::plotDispEsts(dds_wald)
dev.off()

# Median of ratios normalized counts
# Same results if it is Wald or LRT
write_tsv(
  DESeq2::counts(dds_wald, normalized = TRUE),
  outdir_path,
  "median_ratios_normalized_counts.tsv"
)

# PCA with variance-stabilised counts
# Same results if it is Wald or LRT
message("Plotting PCA")
vsd <- DESeq2::vst(dds_wald, blind = FALSE)
write_tsv(SummarizedExperiment::assay(vsd), outdir_path, "vst_counts.tsv")

for (ntop in c(100, 1000, 10000, length(vsd))) {
  if (ntop >= length(vsd)) {
    ntop <- length(vsd)
  }
  plot_pca(vsd, ntop)
  pca_file_path <- paste0(outdir_path, "/pca_ntop_", ntop, ".png")
  ggsave(pca_file_path, dpi = 600)
}

# Heatmap of the sample-to-sample distances
message("Plotting samples heatmap")
png(file.path(outdir_path, "heatmap_samples.png"),
  width = 8, height = 7, units = "in", res = 600
)
plot_heatmap_samples(vsd)
dev.off()
