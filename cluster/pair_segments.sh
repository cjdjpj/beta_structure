#!/bin/bash
#SBATCH --job-name=pair_segments
#SBATCH --time=1-00:00:00                    
#SBATCH --partition=1day
#SBATCH --ntasks=1
#SBATCH --mem=150G

srun python pair_segments.py \
    --input "runs_structured/151"\
    --peak "first_peak"\
    --num_pairs 1\

