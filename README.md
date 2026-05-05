# RNA-seq Analysis — *Lucilia sericata* Thermal Stress

Snakemake pipeline for the L4 larval thermal stress RNA-seq experiment.  
Tests whether genetic diversity (HD vs LD) modulates transcriptomic plasticity under acute cold (4°C) and heat (38°C) stress.

**48 individual samples** | 6 lines | 3 temperatures | 150 bp paired-end | strand-specific (dUTP)

---

## Quick start on NeSI

```bash
# 1. Clone the repository
git clone <repo-url>
cd rna-seq-analysis

# 2. Place raw FASTQ files
#    Naming must match sample_sheet.csv: {sample_id}_R{1,2}.fastq.gz
cp /path/to/sequencing/*.fastq.gz data/raw/

# 3. Add reference files (see Resources section below)

# 4. Update project code and email in scripts/submit_pipeline.sl

# 5. Create the conda environment
module load Miniconda3
conda env create -f envs/rnaseq.yaml
conda activate rnaseq_env

# 6. Dry-run to verify the DAG
snakemake -n --configfile config/config.yaml

# 7. Submit the pipeline
sbatch scripts/submit_pipeline.sl
```

---

## Directory structure

```
.
├── config/config.yaml          pipeline parameters
├── data/
│   ├── sample_sheet.csv        sample metadata (48 samples)
│   └── raw/                    FASTQ files (not tracked by git)
├── docs/sample_tree.txt        full sample/output path reference
├── envs/rnaseq.yaml            conda environment
├── profiles/nesi/config.yaml   Snakemake SLURM profile
├── resources/                  genome, index, adapters (see below)
├── scripts/
│   ├── submit_pipeline.sl      SLURM submission script
│   ├── deseq2_analysis.R
│   ├── go_enrichment.R
│   └── wgcna_analysis.R
└── workflow/
    ├── Snakefile
    └── rules/                  qc | trim | align | count | deseq2
```

---

## Resources (not tracked by git)

Download and place in `resources/` before running:

| File | Source | Path |
|------|--------|------|
| *L. sericata* genome | NCBI/local | `resources/genome/Lucilia_sericata.fna` |
| GTF annotation | NCBI/local | `resources/genome/Lucilia_sericata.gtf` |
| TruSeq adapters | Trimmomatic bundled | `resources/adapters/TruSeq3-PE-2.fa` |
| Gene→GO mapping | Functional annotation | `resources/gene2go.map` |

The HISAT2 index (`resources/hisat2_index/`) is built automatically by the pipeline on first run.

---

## Pipeline steps

| Step | Tool | Key parameters |
|------|------|----------------|
| Raw QC | FastQC 0.11.9 | — |
| Aggregated QC | MultiQC 1.14 | — |
| Adapter trimming | Trimmomatic 0.39 | ILLUMINACLIP TruSeq3-PE-2, MINLEN 36 |
| Alignment | HISAT2 2.2.1 | `--rna-strandness RF --dta` |
| Sort/index | SAMtools 1.17 | — |
| Read counting | featureCounts 2.0.3 | `-s 2 -p -B -C` (reverse-stranded, PE) |
| Differential expression | DESeq2 1.40 | FDR < 0.05, \|LFC\| ≥ 1.0 |
| GO enrichment | topGO 2.52 | Fisher's exact, weight01 algorithm |
| Co-expression network | WGCNA 1.72 | signed hybrid, auto soft-power |

---

## Experimental design

| Line | Diversity | Type | n per temperature |
|------|-----------|------|-------------------|
| LD9-A | LD | Pseudo-replicate (G9 inbred) | 3 |
| LD9-B | LD | Pseudo-replicate (G9 inbred) | 3 |
| LD9-C | LD | Pseudo-replicate (G9 inbred) | 2 |
| FC-HD1 | HD | Field-collected (independent) | 2 |
| FC-HD2 | HD | Field-collected (independent) | 3 |
| FC-HD3 | HD | Field-collected (independent) | 3 |

Temperatures: control 22°C | cold 4°C | heat 38°C | 3-hour exposure

DESeq2 model: `~ line + temperature + diversity:temperature`

**Design note:** LD sub-cages are pseudo-replicates from a single G9 source and are treated as cage-level technical replicates in the model — not independent diversity replicates. See `docs/sample_tree.txt` for full discussion.

---

## Contrasts tested

- `cold_vs_control` — transcriptomic cold stress response
- `heat_vs_control` — transcriptomic heat stress response
- `heat_vs_cold` — asymmetry between stress directions
- `HD_vs_LD_control` — baseline diversity effect
- `HD_vs_LD_cold` — diversity × cold stress interaction
- `HD_vs_LD_heat` — diversity × heat stress interaction

---

## Citation / Methods reference

See the methods section document (`RNA_seq_Methods_v3_Revised.docx`) for full experimental details, QC thresholds (RIN ≥ 7), and statistical rationale.
