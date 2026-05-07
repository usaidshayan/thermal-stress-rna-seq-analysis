#!/usr/bin/env Rscript
# GO enrichment analysis using topGO
# runs on each DESeq2 contrast separately - called once per contrast by snakemake

suppressPackageStartupMessages({
  library(topGO)
  library(ggplot2)
  library(optparse)
})

option_list <- list(
  make_option("--deseq2_results", type = "character"),
  make_option("--contrast",       type = "character"),
  make_option("--outdir",         type = "character"),
  make_option("--fdr",            type = "double", default = 0.05),
  make_option("--go_fdr",         type = "double", default = 0.05,
              help = "FDR threshold for GO term significance [default %default]"),
  make_option("--gene2go",        type = "character", default = "resources/gene2go.map",
              help = "Two-column TSV: gene_id<tab>GO:XXXXXXX (one row per annotation)")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# load DESeq2 results for this contrast and set up gene universe
# drop NAs - genes that DESeq2 couldn't test (outliers, low counts)
res <- read.csv(opt$deseq2_results, stringsAsFactors = FALSE)
res <- res[!is.na(res$padj), ]

# topGO wants a named factor: 1 = significant, 0 = not
all_genes  <- setNames(res$padj, res$gene_id)
sig_genes  <- factor(as.integer(all_genes < opt$fdr))
names(sig_genes) <- res$gene_id

# gene2go mapping file - needs to be generated from the L. sericata annotation
# format is: gene_id\tGO:0000001,GO:0000002 etc.
# if this file is missing the script will fail here with a helpful message
if (!file.exists(opt$gene2go)) {
  stop("gene2go mapping file not found: ", opt$gene2go,
       "\nGenerate it from the L. sericata annotation (see resources/README.md)")
}
gene2go <- readMappings(opt$gene2go)

# runs topGO for a given ontology namespace (BP, MF, or CC)
# using weight01 algorithm which accounts for GO graph structure - reduces redundancy vs classic fisher
run_topgo <- function(ontology) {
  go_data <- new("topGOdata",
                 ontology    = ontology,
                 allGenes    = sig_genes,
                 annot       = annFUN.gene2GO,
                 gene2GO     = gene2go,
                 nodeSize    = 10)  # min 10 genes per term to reduce noise from rare annotations

  result_fisher <- runTest(go_data, algorithm = "weight01", statistic = "fisher")

  n_sig <- sum(score(result_fisher) < opt$go_fdr)
  if (n_sig == 0) {
    message(sprintf("  No significant GO terms (%s) at FDR < %.2f", ontology, opt$go_fdr))
    return(NULL)
  }

  tab <- GenTable(go_data, fisher = result_fisher,
                  orderBy = "fisher", topNodes = min(100, n_sig + 20))
  tab$ontology  <- ontology
  tab$contrast  <- opt$contrast
  tab$fisher_p  <- as.numeric(tab$fisher)
  tab
}

# run for all three ontologies and combine
results_list <- lapply(c("BP", "MF", "CC"), run_topgo)
results_all  <- do.call(rbind, Filter(Negate(is.null), results_list))

# write empty files if nothing significant - snakemake needs the outputs to exist either way
if (is.null(results_all) || nrow(results_all) == 0) {
  message("No significant GO terms for contrast: ", opt$contrast)
  write.csv(data.frame(), file.path(opt$outdir, "go_results.csv"), row.names = FALSE)
  pdf(file.path(opt$outdir, "go_barplot.pdf")); dev.off()
  quit(save = "no")
}

write.csv(results_all, file.path(opt$outdir, "go_results.csv"), row.names = FALSE) #output the results to a csv file

# bar plot of top 20 BP terms - just biological process for now, MF/CC are in the CSV
bp_sig <- subset(results_all, ontology == "BP" & fisher_p < opt$go_fdr)
bp_sig <- bp_sig[order(bp_sig$fisher_p), ][seq_len(min(20, nrow(bp_sig))), ]
bp_sig$Term <- factor(bp_sig$Term, levels = rev(bp_sig$Term))

if (nrow(bp_sig) > 0) {
  pdf(file.path(opt$outdir, "go_barplot.pdf"), width = 10, height = 7)
  print(
    ggplot(bp_sig, aes(x = Term, y = -log10(fisher_p), fill = fisher_p)) +
      geom_col() +
      coord_flip() +
      scale_fill_viridis_c(option = "plasma", direction = -1, name = "p-value") +
      labs(
        title    = sprintf("GO Biological Process enrichment\n%s", opt$contrast),
        x        = NULL,
        y        = expression(-log[10](p-value))
      ) +
      theme_bw(base_size = 11)
  )
  dev.off()
} else {
  pdf(file.path(opt$outdir, "go_barplot.pdf")); dev.off()
}

message("GO enrichment complete: ", opt$contrast)
