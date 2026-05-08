#!/usr/bin/env Rscript
# WGCNA - weighted gene co-expression network analysis
# looking for gene modules that correlate with thermal treatment and/or diversity

suppressPackageStartupMessages({
  library(WGCNA)
  library(ggplot2)
  library(pheatmap)
  library(optparse)
})

option_list <- list(
  make_option("--vst_counts",       type = "character"),
  make_option("--sample_sheet",     type = "character"),
  make_option("--outdir",           type = "character", default = "results/wgcna/"),
  make_option("--soft_power",       type = "character", default = "auto"),
  make_option("--min_module_size",  type = "integer",   default = 30),
  make_option("--merge_cut_height", type = "double",    default = 0.25),
  make_option("--threads",          type = "integer",   default = 8)
)
opt <- parse_args(OptionParser(option_list = option_list))

options(stringsAsFactors = FALSE)
enableWGCNAThreads(opt$threads)
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# load VST counts from DESeq2 - blind=FALSE VST is appropriate here since we're not doing QC
vst <- read.csv(opt$vst_counts, row.names = 1, check.names = FALSE)
meta <- read.csv(opt$sample_sheet, stringsAsFactors = FALSE)
meta <- meta[match(colnames(vst), meta$sample_id), ]

# WGCNA wants samples as rows, genes as columns - opposite of how the count matrix comes
expr_data <- t(as.matrix(vst))

# drop bottom 25% of genes by MAD - low variance genes don't contribute useful network structure
# and keeping them just slows everything down
gene_mad   <- apply(expr_data, 2, mad)
keep_genes <- gene_mad >= quantile(gene_mad, 0.25)
expr_data  <- expr_data[, keep_genes]
message(sprintf("Using %d genes (top 75%% by MAD) for WGCNA", ncol(expr_data)))

# check for bad samples/genes before building the network
gsg <- goodSamplesGenes(expr_data, verbose = 0)
if (!gsg$allOK) {
  expr_data <- expr_data[gsg$goodSamples, gsg$goodGenes]
  message("Removed outlier samples/genes flagged by goodSamplesGenes()")
}

# --- soft-thresholding power ---
# this is the key WGCNA tuning parameter - we want the lowest power that gives R^2 > 0.85
# script will try to pick it automatically, falling back to 12 if it can't
powers <- c(1:20)
sft    <- pickSoftThreshold(expr_data, powerVector = powers, verbose = 0,
                             networkType = "signed hybrid")

if (opt$soft_power == "auto") {
  soft_power <- sft$powerEstimate
  if (is.na(soft_power)) {
    soft_power <- 12   # 12 is a reasonable fallback for signed hybrid networks
    message("Automatic soft power selection failed; using fallback = 12")
  }
} else {
  soft_power <- as.integer(opt$soft_power)
}
message(sprintf("Using soft-thresholding power: %d", soft_power))

# save the soft power plots so we can check the choice looks reasonable
pdf(file.path(opt$outdir, "soft_power_selection.pdf"), width = 9, height = 5)
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft threshold (power)", ylab = "Scale-free R²", type = "n",
     main = "Scale-free fit")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = "red")
abline(h = 0.85, col = "darkgreen", lty = 2)  # 0.85 is the recommended threshold
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft threshold (power)", ylab = "Mean connectivity", type = "n",
     main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")
dev.off()

# blockwiseModules is more memory efficient than the one-shot approach for larger gene sets
# saveTOMs=FALSE because the TOM files get huge and we don't need them downstream
net <- blockwiseModules(
  expr_data,
  power             = soft_power,
  networkType       = "signed hybrid",
  TOMType           = "signed",
  minModuleSize     = opt$min_module_size,
  mergeCutHeight    = opt$merge_cut_height,  # merge modules with similar eigengenes
  numericLabels     = FALSE,
  pamRespectsDendro = FALSE,
  saveTOMs          = FALSE,
  verbose           = 0
)

module_colours <- net$colors
module_eigengenes <- net$MEs

message(sprintf("Detected %d modules (excluding grey/unassigned)",
                length(unique(module_colours)) - 1))

# grey module = unassigned genes, that's fine and expected
pdf(file.path(opt$outdir, "dendrogram_and_colours.pdf"), width = 14, height = 6)
plotDendroAndColors(
  net$dendrograms[[1]],
  module_colours[net$blockGenes[[1]]],
  "Module colours",
  dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05,
  main = "Gene dendrogram and module colours"
)
dev.off()

# --- module-trait correlations ---

# encode traits as binary columns - WGCNA correlation works on numeric data
trait_data <- data.frame(
  cold    = as.integer(meta$temp_label == "cold"),
  heat    = as.integer(meta$temp_label == "heat"),
  control = as.integer(meta$temp_label == "control"),
  HD      = as.integer(meta$diversity == "HD"),
  row.names = meta$sample_id
)

# remove grey module eigengene before correlating - it's just noise
me_clean  <- orderMEs(removeGreyME(module_eigengenes))
mod_trait_cor  <- cor(me_clean, trait_data, use = "p")
mod_trait_pval <- corPvalueStudent(mod_trait_cor, nrow(expr_data))

cor_df <- as.data.frame(mod_trait_cor)
cor_df$module <- rownames(cor_df)
write.csv(cor_df, file.path(opt$outdir, "module_trait_correlations.csv"), row.names = FALSE)

gene_df <- data.frame(
  gene_id = colnames(expr_data),
  module  = module_colours
)
write.csv(gene_df, file.path(opt$outdir, "module_gene_assignments.csv"), row.names = FALSE)

# heatmap showing r and p-value for each module x trait combination
text_matrix <- paste0(
  signif(mod_trait_cor, 2), "\n(",
  signif(mod_trait_pval, 1), ")"
)
dim(text_matrix) <- dim(mod_trait_cor)

pdf(file.path(opt$outdir, "module_trait_heatmap.pdf"), width = 8, height = max(6, nrow(me_clean) * 0.4))
par(mar = c(5, 10, 3, 3))
labeledHeatmap(
  Matrix     = mod_trait_cor,
  xLabels    = colnames(trait_data),
  yLabels    = rownames(mod_trait_cor),
  ySymbols   = rownames(mod_trait_cor),
  colorLabels = FALSE,
  colors     = blueWhiteRed(50),
  textMatrix = text_matrix,
  setStdMargins = FALSE,
  cex.text   = 0.7,
  zlim       = c(-1, 1),
  main       = "Module–trait correlations\n(Pearson r, p-value)"
)
dev.off()

message("WGCNA analysis complete.")
