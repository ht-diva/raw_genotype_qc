#!/usr/bin/env bash

#SBATCH --job-name=raw_genotype_qc
#SBATCH --output=%x_%j.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=username@fht.org
#SBATCH --partition=cpuq
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=08:00:00

set -eo pipefail

# Move to the directory from which sbatch was executed.
cd "${SLURM_SUBMIT_DIR}"

# Configure Singularity/Apptainer for the HDS cluster.
source /exchange/healthds/singularity_functions

# Activate Snakemake.
source ~/.bashrc
conda activate /exchange/healthds/software/envs/snakemake

echo "Pipeline directory: $(pwd)"
echo "Execution host: $(hostname)"
echo "Start time: $(date)"

make run

echo "Completion time: $(date)"