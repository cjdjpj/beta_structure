#!/bin/bash

INPUT="runs/r001"

sbatch dist.sh "$INPUT"
sbatch frac_iden_blk.sh "$INPUT"
sbatch frac_clonal.sh "$INPUT"
sbatch rd.sh "$INPUT"
