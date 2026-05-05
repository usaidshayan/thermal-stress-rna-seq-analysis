#!/usr/bin/env Rscript
# DESeq2 differential expression - L. sericata thermal stress
# HD vs LD comparisons + temperature contrasts

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(optparse)
})

option_list <- list(
  make_option("--counts",       type = "character"),
  make_option("--sample_sheet", type = "character"),
  make_option("--outdir",       type = "character", default = "results/deseq2/"),
  make_option("--fdr",          type = "double",    default = 0.05),
  make_option("--lfc",          type = "double",    default = 1.0),
  make_option("--threads",      type = "integer",   default = 4)
)
opt <- parse_args(OptionParser(option_list = option_list))

BiocParallel::register(BiocParallel::MulticoreParam(opt$threads))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# --- load data ---

raw_counts <- read.delim(opt$counts, comment.char = "#", check.names = FALSE)

# featureCounts puts metadata in cols 1-6, actual counts start at col 7
count_matrix <- as.matrix(raw_counts[, 7:ncol(raw_counts)])
rownames(count_matrix) <- raw_counts$Geneid

# column names come out as full BAM paths, strip them back to sample IDs
colnames(count_matrix) <- gsub(".*/|.sorted.bam", "", colnames(count_matrix))

meta <- read.csv(opt$sample_sheet, stringsAsFactors = FALSE) #feeding the sample sheet 
meta <- meta[match(colnames(count_matrix), meta$sample_id), ]
stopifnot(all(meta$sample_id == colnames(count_matrix)))  # paranoia check - order must match

# set reference levels: control is baseline for temperature, LD is baseline for diversity
meta$temperature <- factor(meta$temp_label, levels = c("control", "cold", "heat"))
meta$diversity    <- factor(meta$diversity,  levels = c("LD", "HD"))
meta$line         <- factor(meta$line)

# model: line as a blocking factor to absorb cage/line variance
# then temperature main effect + diversity:temperature interaction
# the interaction is what we actually care about - does diversity change the thermal response?
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData   = meta,
  design    = ~ line + temperature + diversity:temperature
)

# pre-filter genes with very low counts - speeds things up and reduces multiple testing burden
# keeping genes with >=10 counts in at least n_min samples (n_min = smallest treatment group)
n_min <- min(table(meta$temp_label))
keep  <- rowSums(counts(dds) >= 10) >= n_min
dds   <- dds[keep, ]
message(sprintf("Retained %d genes after low-count filtering", sum(keep)))

dds <- DESeq(dds, parallel = TRUE)

# save normalised counts and VST - VST is used downstream for PCA and WGCNA
norm_counts <- counts(dds, normalized = TRUE)
write.csv(as.data.frame(norm_counts),
          file.path(opt$outdir, "normalised_counts.csv"))

vst_counts  <- assay(vst(dds, blind = FALSE))
write.csv(as.data.frame(vst_counts),
          file.path(opt$outdir, "vst_counts.csv"))

# pca 

pca_data <- plotPCA(vst(dds, blind = FALSE),
                    intgroup = c("diversity", "temperature"),
                    returnData = TRUE)
pct_var <- round(100 * attr(pca_data, "percentVar"))

pdf(file.path(opt$outdir, "pca_plot.pdf"), width = 8, height = 6)
ggplot(pca_data, aes(PC1, PC2, colour = temperature, shape = diversity, label = name)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(size = 2.5, max.overlaps = 20) +
  scale_colour_manual(values = c(control = "#2196F3", cold = "#9C27B0", heat = "#F44336")) +
  labs(
    x = sprintf("PC1: %d%% variance", pct_var[1]),
    y = sprintf("PC2: %d%% variance", pct_var[2]),
    title = "PCA — VST-normalised gene expression",
    subtitle = "L. sericata L4 thermal stress (n = 48)"
  ) +
  theme_bw(base_size = 12)
dev.off()

# sample distance heatmap - good sanity check for clustering by treatment
samp_dists <- dist(t(vst_counts))
samp_dist_mat <- as.matrix(samp_dists)
ann_col <- data.frame(
  Temperature = meta$temperature,
  Diversity   = meta$diversity,
  Line        = meta$line,
  row.names   = meta$sample_id
)

pdf(file.path(opt$outdir, "sample_distance_heatmap.pdf"), width = 12, height = 10)
pheatmap(samp_dist_mat,
         clustering_distance_rows = samp_dists,
         clustering_distance_cols = samp_dists,
         annotation_col = ann_col,
         color = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         main = "Sample-to-sample Euclidean distances (VST)")
dev.off()

# --- contrasts ---

# helper function - writes results table, significant gene list, and MA plot for a given contrast
save_results <- function(res, contrast_name) {
  out_dir <- file.path(opt$outdir, contrast_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  res_df <- as.data.frame(res) |>
    tibble::rownames_to_column("gene_id") |>
    dplyr::arrange(padj)

  write.csv(res_df, file.path(out_dir, "results_table.csv"), row.names = FALSE)

  # significant = FDR threshold AND absolute LFC threshold
  sig <- subset(res_df, padj < opt$fdr & abs(log2FoldChange) >= opt$lfc)
  write.csv(sig, file.path(out_dir, "significant_genes.csv"), row.names = FALSE)
  message(sprintf("%s: %d DEGs (FDR < %.2f, |LFC| >= %.1f)",
                  contrast_name, nrow(sig), opt$fdr, opt$lfc))

  pdf(file.path(out_dir, "ma_plot.pdf"), width = 7, height = 5)
  DESeq2::plotMA(res, ylim = c(-8, 8), main = contrast_name, alpha = opt$fdr)
  dev.off()
}

# temperature contrasts first - these are straightforward main effects
# LD sub-cages are treated as blocked replicates, which the model handles via the line term
save_results(
  results(dds, contrast = c("temperature", "cold", "control"), alpha = opt$fdr),
  "cold_vs_control"
)
save_results(
  results(dds, contrast = c("temperature", "heat", "control"), alpha = opt$fdr),
  "heat_vs_control"
)
save_results(
  results(dds, contrast = c("temperature", "heat", "cold"), alpha = opt$fdr),
  "heat_vs_cold"
)

# now the diversity x temperature interaction contrasts
# these pull from the interaction terms in the model - check resultsNames(dds) if anything's missing
for (temp in c("control", "cold", "heat")) {
  coef_name <- paste0("diversityHD.temperature", temp)
  if (coef_name %in% resultsNames(dds)) {
    save_results(
      results(dds, name = coef_name, alpha = opt$fdr),
      paste0("HD_vs_LD_", temp)
    )
  } else {
    message("Contrast term not found: ", coef_name, " — skipping")
  }
}

message("DESeq2 analysis complete.")
