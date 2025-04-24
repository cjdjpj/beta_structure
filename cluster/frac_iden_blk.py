import tskit
import pickle
import argparse
import numpy as np
from itertools import combinations

parser = argparse.ArgumentParser(
                    prog='frac_iden_blk')
parser.add_argument('--input', type=str, default="output")
parser.add_argument('--blk_size', type=int, default=10)

args = parser.parse_args()

# open tree sequence
mts = tskit.load(args.input)

pairs = list(combinations(range(mts.num_samples), 2))

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
for (i,j) in pairs:
    frac_iden_blk.append(prop_identical_blk(genotypes[:, i], genotypes[:, j]))

with open(args.input + "_frac_iden_blk", 'wb') as file:
    pickle.dump(frac_iden_blk, file)

print("---fraction of identical blocks computed---")
