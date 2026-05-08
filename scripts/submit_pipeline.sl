#!/bin/bash
# this script submits with: sbatch scripts/submit_pipeline.sl 
#
# barely uses any mem/cpu, job is just the snakemake orchestrator
# actual compute occurs in sub-jobs that snakemake submits to slurm autonomously

#SBATCH --job-name=rnaseq_lsericata
#SBATCH --account=YOUR_PROJECT_CODE          # <-- you'd have to update this to your project code 
#SBATCH --time=72:00:00                      # 72h should be sufficient for the entire pipeline
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2
#SBATCH --partition=large
#SBATCH --output=logs/slurm/orchestrator_%j.out
#SBATCH --error=logs/slurm/orchestrator_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=YOUR_EMAIL@example.com   # <-- update this as well

module purge
module load Miniconda3
source activate rnaseq_env   # build this first with: conda env create -f envs/rnaseq.yaml, this caused me a lot of trouble

mkdir -p logs/slurm

# snakemake reads the profile for SLURM settings and submits each rule as its own job
# --jobs 50 = max 50 jobs running at once (adjust if your project has tighter limits)
# --rerun-incomplete picks up where we left off if the job gets killed 
# --keep-going essentially means if one sample fails the rest still run
snakemake \
    --profile profiles/nesi \
    --configfile config/config.yaml \
    --jobs 50 \
    --rerun-incomplete \
    --keep-going \
    --latency-wait 60 \     # NeSI lustre filesystem can be slow to register new files, things i find out everyday
    --printshellcmds \
    2>&1 | tee logs/slurm/snakemake_$(date +%Y%m%d_%H%M%S).log