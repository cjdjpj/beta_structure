#!/bin/bash
#SBATCH --job-name=n_snp
#SBATCH --time=1-00:00:00                    
#SBATCH --partition=1day
#SBATCH --ntasks=1
#SBATCH --mem=150G
#SBATCH --array=0-8

n_values=(2 3 4 5 6 7 8 9 11)

n=${n_values[$SLURM_ARRAY_TASK_ID]}

srun python n_snp.py \
    --input "runs/r004"\
    --n "$n"

