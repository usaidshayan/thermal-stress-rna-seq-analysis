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