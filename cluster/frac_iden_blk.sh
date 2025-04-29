#!/bin/bash
#SBATCH --job-name=frac_iden_blk
#SBATCH --time=1-00:00:00                    
#SBATCH --partition=1day
#SBATCH --ntasks=1
#SBATCH --mem=150G

srun python frac_iden_blk.py\
    --input "$1"\
    --blk_size 500
