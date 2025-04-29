#!/bin/bash
#SBATCH --job-name=frac_clonal
#SBATCH --time=3-00:00:00                    
#SBATCH --partition=3day
#SBATCH --ntasks=1
#SBATCH --mem=150G

srun python frac_clonal.py\
    --input "$1"\
