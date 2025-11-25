#!/bin/bash
#SBATCH --job-name=biggest_burst
#SBATCH --time=3-00:00:00                    
#SBATCH --partition=3day
#SBATCH --ntasks=1
#SBATCH --mem=150G

srun python biggest_burst.py\
    --input "$1"
