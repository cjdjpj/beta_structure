#!/bin/bash

INPUT="r001"

sbatch dist.sh "runs/$INPUT"
sbatch frac_iden_blk.sh "runs/$INPUT"
sbatch frac_clonal.sh "runs/$INPUT"
