#!/bin/bash
#SBATCH --job-name=sim
#SBATCH --partition=4hrs
#SBATCH --mem=4G
#SBATCH --ntasks=1

srun python sim.py\
    --output "runs/r001"\
    --length 5000000\
    --Ne 25000000\
    --track_length 5000\
    --nsample 500\
    --mu 0.0000000006\
    --r_m 0.0\
    --model "kingman"\
