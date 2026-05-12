# QC - fastqc on raw reads then again after trimming, multiqc to aggregate
# NeSI modules: FastQC/0.11.9, MultiQC/1.14

rule fastqc_raw:
    input:
        r1 = get_r1,
        r2 = get_r2,
    output:
        r1_html = "results/qc/fastqc_raw/{sample}_R1_fastqc.html",
        r2_html = "results/qc/fastqc_raw/{sample}_R2_fastqc.html",
        r1_zip  = "results/qc/fastqc_raw/{sample}_R1_fastqc.zip",
        r2_zip  = "results/qc/fastqc_raw/{sample}_R2_fastqc.zip",
    log:
        "logs/fastqc_raw/{sample}.log"
    threads: 2
    resources:
        mem_mb   = 4000, # i've found higher memory usage to not benefit runtime
        runtime  = 30,
        partition = "large",
    shell:
        """
        module load FastQC/0.11.9 # you may want to double check the version here, has a tendancy of breaking
        fastqc --threads {threads} --outdir results/qc/fastqc_raw/ \
            {input.r1} {input.r2} 2> {log}
        # fastqc names outputs after the input file, not the wildcard sample name
        # so we rename them here to match what snakemake expects as outputs
        mv results/qc/fastqc_raw/$(basename {input.r1} .fastq.gz)_fastqc.html {output.r1_html}
        mv results/qc/fastqc_raw/$(basename {input.r2} .fastq.gz)_fastqc.html {output.r2_html}
        mv results/qc/fastqc_raw/$(basename {input.r1} .fastq.gz)_fastqc.zip  {output.r1_zip}
        mv results/qc/fastqc_raw/$(basename {input.r2} .fastq.gz)_fastqc.zip  {output.r2_zip}