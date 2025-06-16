#!/bin/bash
#SBATCH --job-name=rd
#SBATCH --time=1-00:00:00                    
#SBATCH --partition=1day
#SBATCH --ntasks=1
#SBATCH --mem=150G

srun python rd.py\
    --input "runs/r001"
