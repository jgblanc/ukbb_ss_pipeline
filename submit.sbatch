#!/bin/bash
#SBATCH --job-name=ukbb_testrun_big
#SBATCH --output=logs/ukbb_small1.out
#SBATCH --error=logs/ukbb_small1.err
#SBATCH --time=04:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=14
#SBATCH --mem-per-cpu=4000

module load python/cpython-3.7.0
module load R

echo "SLURM_JOBID="$SLURM_JOBID
cat snakefile
snakemake --cores all --rerun-incomplete









