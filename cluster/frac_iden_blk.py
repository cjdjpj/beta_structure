import json
import tskit
import pickle
import random
import argparse
import numpy as np
from itertools import combinations

parser = argparse.ArgumentParser(
                    prog='frac_iden_blk')
parser.add_argument('--input', type=str, default="output")
parser.add_argument('--num_pairs', type=int, default=0)
parser.add_argument('--blk_size', type=int, default=10)

args = parser.parse_args()

# open files
mts = tskit.load(args.input)
with open(args.input + ".json", "r") as file:
    params = json.load(file)

# choose pairs
pairs = list(combinations(range(params["nsample"]), 2))
random.seed(42)
random_pair_indices = random.sample(range(len(pairs)), args.num_pairs)

# compute
sites = [int(site.position) for site in mts.sites()]

# tally number of sites per block
block_indices = np.floor_divide(sites, args.blk_size)
muts_per_blk = np.bincount(block_indices)
muts_per_blk = muts_per_blk[muts_per_blk != 0]

def prop_identical_blk(s1, s2):
    i = 0
    matches = 0
    for muts in muts_per_blk:
        if np.array_equal(s1[i : i + muts], s2[i : i + muts]):
            matches += 1
        i += muts
    return matches / len(muts_per_blk)

frac_iden_blk= []
genotypes = mts.genotype_matrix()
for pair_index in random_pair_indices:
    (i,j) = pairs[pair_index]
    frac_iden_blk.append(prop_identical_blk(genotypes[:, i], genotypes[:, j]))

with open(args.input + "_frac_iden_blk", 'wb') as file:
    pickle.dump(frac_iden_blk, file)

print("---fraction of identical blocks computed---")
