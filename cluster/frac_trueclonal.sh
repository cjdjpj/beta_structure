#!/bin/bash
#SBATCH --job-name=frac_trueclonal
#SBATCH --time=3-00:00:00                    
#SBATCH --partition=3day
#SBATCH --ntasks=1

srun python frac_trueclonal.py\
    --input "$1"\
