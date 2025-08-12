#!/bin/bash
#SBATCH --job-name=n_snp_entropy
#SBATCH --time=1-00:00:00                    
#SBATCH --partition=1day
#SBATCH --ntasks=1
#SBATCH --mem=150G

srun python n_snp_entropy.py \
    --input "runs/r002"\

