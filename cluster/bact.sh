#!/bin/bash
#SBATCH --job-name=bact_coalescent            # Job name
#SBATCH --mem=2G                              # Memory request
#SBATCH --ntasks=1                            # Number of tasks
#SBATCH --cpus-per-task=1                     # Number of CPU cores per task

srun python bact.py\
    --output "runs/r001"\
    --length 5000000\
    --Ne 25000000\
    --track_length 5000\
    --nsample 500\
    --mu 0.0000000006\
    --r_m 0.0\
    --model "kingman"\
