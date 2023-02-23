#!/bin/bash
#SBATCH --job-name=get_snp_num
#SBATCH --output=logs/get_snp_num.out
#SBATCH --error=logs/get_snp_num.err
#SBATCH --time=06:00:00
#SBATCH --partition=build
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=14
#SBATCH --mem-per-cpu=2000

module load python/cpython-3.7.0
module load R

echo "SLURM_JOBID="$SLURM_JOBID
cat snakefile
snakemake --cores all --rerun-incomplete --keep-going









