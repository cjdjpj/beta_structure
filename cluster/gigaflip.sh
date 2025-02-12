#!/bin/bash
#SBATCH --job-name=gigaflip
#SBATCH --partition=4hrs
#SBATCH --mem=4G
#SBATCH --ntasks=1

srun ./gigaflip
