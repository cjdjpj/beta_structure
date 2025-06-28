#!/bin/bash
#SBATCH --job-name=mass_sim
#SBATCH --partition=1day
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=250
#SBATCH --mem=15G
#SBATCH --array=0-39
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

START=$((SLURM_ARRAY_TASK_ID * 250))
END=$((START + 249))

for i in $(seq $START $END); do
    python mass_sim.py \
        --output "mass_sim_results/k_0_${i}" \
        --length 5000000 \
        --track_length 5000 \
        --nsample 100 \
        --mu 0.025 \
        --r_m 0.0 \
        --model "kingman" \
        --pi 0.03 \
    &> /dev/null &
done

wait

