#!/bin/bash
#SBATCH --job-name=fract_ident_blk
#SBATCH --time=1-00:00:00                    
#SBATCH --partition=1day
#SBATCH --ntasks=1

srun python frac_iden_blk.py\
    --input "$1"\
    --num_pairs 28800\
    --blk_size 1000
