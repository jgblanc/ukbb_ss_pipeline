#!/bin/bash
#SBATCH --job-name=ukbb_all_weighted_strict
#SBATCH --output=logs/ukbb_all_weighted_strict.out
#SBATCH --error=logs/ukbb_all_weighted_strict.err
#SBATCH --time=12:00:00
#SBATCH --partition=build
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=14
#SBATCH --mem-per-cpu=4000

module load python/cpython-3.7.0
module load R

echo "SLURM_JOBID="$SLURM_JOBID
cat snakefile
snakemake --cores all --rerun-incomplete --keep-going









