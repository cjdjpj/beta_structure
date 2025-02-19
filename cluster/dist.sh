#!/bin/bash
#SBATCH --job-name=dist
#SBATCH --time=04:00:00                    
#SBATCH --partition=4hrs
#SBATCH --ntasks=1

srun python dist.py\
    --input "$1"
