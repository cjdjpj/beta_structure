#!/bin/bash
#SBATCH --job-name=n_snp
#SBATCH --time=1-00:00:00                    
#SBATCH --partition=1day
#SBATCH --ntasks=1
#SBATCH --mem=150G
#SBATCH --array=0-16

n_values=(1 2 3 4 5 6 7 8 9 11 13 15 17 19 21 23 25)

n=${n_values[$SLURM_ARRAY_TASK_ID]}

srun python n_snp.py \
    --input "runs_inf_sites/is01"\
    --n "$n"

