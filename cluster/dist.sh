#!/bin/bash
#SBATCH --job-name=dist
#SBATCH --time=3-00:00:00                    
#SBATCH --partition=3day
#SBATCH --ntasks=1
#SBATCH --mem=150G

srun python dist.py\
    --input "$1"
