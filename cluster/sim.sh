#!/bin/bash
#SBATCH --job-name=sim
#SBATCH --partition=4hrs
#SBATCH --ntasks=1
#SBATCH --mem=150G

srun python sim.py\
    --output "runs/r001"\
    --length 5000000\
    --track_length 5000\
    --nsample 100\
    --mu 0.015\
    --r 0.0\
    --KT_2 1.0\
    --model "kingman"\
