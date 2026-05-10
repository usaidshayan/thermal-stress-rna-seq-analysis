# downstream stats - DESeq2 -> GO enrichment -> WGCNA
# all running in R, NeSI module: R/4.3.1-gimkl-2022a
# GO enrichment runs once per contrast so snakemake can parallelise the 6 contrasts

CONTRASTS = [
    "cold_vs_control",
    "heat_vs_control",
    "heat_vs_cold",
    "HD_vs_LD_control",  # diversity effect within each temp treatment
    "HD_vs_LD_cold",
    "HD_vs_LD_heat",
]

# DESeq2 produces all contrast tables in one go (can't easily split per-contrast without rerunning the model)
rule deseq2:
    input:
        counts       = "results/counts/all_samples_counts.txt",
        sample_sheet = config["sample_sheet"],
    output:
        expand(
            "results/deseq2/{contrast}/results_table.csv",
            contrast=CONTRASTS
        ),
        "results/deseq2/normalised_counts.csv",
        "results/deseq2/vst_counts.csv",   # VST used for PCA and WGCNA
        "results/deseq2/pca_plot.pdf",
        "results/deseq2/sample_distance_heatmap.pdf",
    log:
        "logs/deseq2/deseq2.log"
    threads: 4
    resources:
        mem_mb   = 32000,
        runtime  = 120,
        partition = "large",
    shell:
        """
        module load R/4.3.1-gimkl-2022a
        Rscript scripts/deseq2_analysis.R \
            --counts {input.counts} \
            --sample_sheet {input.sample_sheet} \
            --outdir results/deseq2/ \
            --fdr {config[deseq2][fdr_threshold]} \
            --lfc {config[deseq2][lfc_threshold]} \
            --threads {threads} \
            2> {log}
        """

# GO enrichment - runs per contrast, snakemake handles the parallelism
rule go_enrichment:
    input:
        deseq2_results = "results/deseq2/{contrast}/results_table.csv",
        sample_sheet   = config["sample_sheet"],
    output:
        "results/go_enrichment/{contrast}/go_results.csv",
        "results/go_enrichment/{contrast}/go_barplot.pdf",
    log:
        "logs/go_enrichment/{contrast}.log"
    threads: 2
    resources:
        mem_mb   = 16000,
        runtime  = 60,
        partition = "large",
    shell:
        """
        module load R/4.3.1-gimkl-2022a
        Rscript scripts/go_enrichment.R \
            --deseq2_results {input.deseq2_results} \
            --contrast {wildcards.contrast} \
            --outdir results/go_enrichment/{wildcards.contrast}/ \
            --fdr {config[deseq2][fdr_threshold]} \
            2> {log}
        """

# WGCNA needs the VST counts from DESeq2 - give it plenty of memory, it can be greedy
rule wgcna:
    input:
        vst_counts   = "results/deseq2/vst_counts.csv",
        sample_sheet = config["sample_sheet"],
    output:
        "results/wgcna/module_trait_correlations.csv",
        "results/wgcna/module_gene_assignments.csv",
        "results/wgcna/dendrogram_and_colours.pdf",
        "results/wgcna/module_trait_heatmap.pdf",
    log:
        "logs/wgcna/wgcna.log"
    params:
        soft_power       = config["wgcna"]["soft_power"],
        min_module_size  = config["wgcna"]["min_module_size"],
        merge_cut_height = config["wgcna"]["merge_cut_height"],
    threads: 8
    resources:
        mem_mb   = 64000,   # well, WGCNA can use a lot of RAM with large gene sets
        runtime  = 240,
        partition = "large",
    shell:
        """
        module load R/4.3.1-gimkl-2022a
        Rscript scripts/wgcna_analysis.R \
            --vst_counts {input.vst_counts} \
            --sample_sheet {input.sample_sheet} \
            --outdir results/wgcna/ \
            --soft_power {params.soft_power} \
            --min_module_size {params.min_module_size} \
            --merge_cut_height {params.merge_cut_height} \
            --threads {threads} \
            2> {log}
        """
