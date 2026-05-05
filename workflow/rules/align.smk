# alignment with HISAT2, then sort+index with samtools
# NeSI modules: HISAT2/2.2.1-gimkl-2022a, SAMtools/1.16.1-GCC-11.3.0
#
# important flags to know about:
#   --dta          tells hisat2 we're going into featureCounts downstream, improves junction handling
#   --rna-strandness RF    library is dUTP reverse-stranded (poly-A selected), RF is correct for this

rule hisat2_index:
    """build the genome index - only runs once, touch file used to track completion"""
    input:
        fasta = config["genome_fasta"],
        gtf   = config["genome_gtf"],
    output:
        touch("resources/hisat2_index/.done")
    log:
        "logs/hisat2_index.log"
    params:
        prefix = config["hisat2_index"],
    threads: 8
    resources:
        mem_mb   = 32000,
        runtime  = 120,
        partition = "large",
    shell:
        """
        module load HISAT2/2.2.1-gimkl-2022a
        # pull splice sites and exons from the GTF first - makes alignment splice-aware
        hisat2_extract_splice_sites.py {input.gtf} > resources/hisat2_index/splice_sites.txt 2>> {log}
        hisat2_extract_exons.py        {input.gtf} > resources/hisat2_index/exons.txt        2>> {log}
        hisat2-build \
            --ss resources/hisat2_index/splice_sites.txt \
            --exon resources/hisat2_index/exons.txt \
            -p {threads} \
            {input.fasta} {params.prefix} 2>> {log}
        """

rule hisat2_align:
    input:
        r1    = "results/trimmed/{sample}_R1_trimmed.fastq.gz",
        r2    = "results/trimmed/{sample}_R2_trimmed.fastq.gz",
        index = "resources/hisat2_index/.done",  # just needs the index to exist
    output:
        bam   = "results/aligned/{sample}.sorted.bam",
        bai   = "results/aligned/{sample}.sorted.bam.bai",
    log:
        hisat2  = "logs/align/{sample}.hisat2.log",
        samtools = "logs/align/{sample}.samtools.log",
    params:
        index      = config["hisat2_index"],
        strandness = config["hisat2"]["strandness"],
        extra      = config["hisat2"]["extra"],
        # read group tags - useful if we ever need to go back and check which sample is which in the BAM
        rg_id      = lambda wildcards: wildcards.sample,
        rg_sm      = lambda wildcards: samples.loc[wildcards.sample, "line"],
        rg_lb      = lambda wildcards: samples.loc[wildcards.sample, "diversity"],
    threads: 8
    resources:
        mem_mb   = 24000,
        runtime  = 120,
        partition = "large",
    shell:
        """
        module load HISAT2/2.2.1-gimkl-2022a SAMtools/1.16.1-GCC-11.3.0
        hisat2 \
            -x {params.index} \
            -1 {input.r1} -2 {input.r2} \
            --rna-strandness {params.strandness} \
            {params.extra} \
            --rg-id {params.rg_id} \
            --rg SM:{params.rg_sm} \
            --rg LB:{params.rg_lb} \
            --rg PL:ILLUMINA \
            -p {threads} \
            --summary-file {log.hisat2} \
            2>> {log.hisat2} \
        | samtools sort -@ {threads} -o {output.bam} - 2> {log.samtools}
        samtools index {output.bam} 2>> {log.samtools}
        """

# flagstat gives mapping rates etc - multiqc picks this up for the alignment QC report
rule samtools_flagstat:
    input:
        bam = "results/aligned/{sample}.sorted.bam",
    output:
        "results/aligned/{sample}.flagstat",
    log:
        "logs/align/{sample}.flagstat.log"
    threads: 2
    resources:
        mem_mb   = 4000,
        runtime  = 20,
        partition = "large",
    shell:
        """
        module load SAMtools/1.16.1-GCC-11.3.0
        samtools flagstat -@ {threads} {input.bam} > {output} 2> {log}
        """
