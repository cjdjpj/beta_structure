#!/bin/bash
#SBATCH --job-name=coal_times
#SBATCH --time=3-00:00:00                    
#SBATCH --partition=3day
#SBATCH --ntasks=1
#SBATCH --mem=150G

srun python coal_times.py\
    --input "$1"
