#!/bin/bash
#SBATCH --job-name=mass_sim_outlier
#SBATCH --partition=1day
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem-per-cpu=35G
#SBATCH --array=0-3

START=$((SLURM_ARRAY_TASK_ID * 50))
END=$((START + 49))

for i in $(seq $START $END); do
    python mass_sim_outlier.py \
        --output "mass_sim_results/b1.1_0.3_${i}" \
        --length 5000000 \
        --track_length 5000 \
        --nsample 100 \
        --mu 0.025 \
        --r_m 0.3 \
        --model "beta" \
        --alpha 1.1 \
        --pi 0.03 \
    &
done

wait
