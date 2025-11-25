#!/bin/bash
#SBATCH --job-name=mass_sim_arity
#SBATCH --partition=1day
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem-per-cpu=35G
#SBATCH --array=0-3

START=$((SLURM_ARRAY_TASK_ID * 50))
END=$((START + 49))

for i in $(seq $START $END); do
    python mass_sim_arity.py \
        --output "mass_sim_results/b1.1_45_${i}" \
        --length 5000000 \
        --track_length 5000 \
        --nsample 100 \
        --mu 0.015 \
        --r 0.0045\
        --KT_2 1.0\
        --model "beta" \
        --alpha 1.1 \
    &
done

wait
