# featureCounts - counts reads per gene across all samples at once
# much easier to get a single matrix than merging 48 individual count files later
# NeSI module: Subread/2.0.3-GCC-11.3.0
#
# key flags (these matter for strand-specific PE data):
#   -s 2   reverse stranded - must match library prep (dUTP method = RF = -s 2)
#   -p     paired end mode
#   -B     only count fragments where both reads map - reduces noise
#   -C     discard chimeric pairs (shouldn't be many but good practice)

rule featurecounts:
    """all 48 BAMs go in together so we get one count matrix out"""
    input:
        bams = expand("results/aligned/{sample}.sorted.bam", sample=SAMPLES),
        gtf  = config["genome_gtf"],
    output:
        counts  = "results/counts/all_samples_counts.txt",
        summary = "results/counts/all_samples_counts.txt.summary",  # per-sample assignment stats
    log:
        "logs/featurecounts/featurecounts.log"
    params:
        extra = config["featurecounts"]["extra"],
    threads: 8
    resources:
        mem_mb   = 16000,
        runtime  = 60,
        partition = "large",
    shell:
        """
        module load Subread/2.0.3-GCC-11.3.0
        featureCounts \
            {params.extra} \
            -a {input.gtf} \
            -o {output.counts} \
            {input.bams} \
            2> {log}
        """
