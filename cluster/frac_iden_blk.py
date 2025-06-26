import tskit
import pickle
import argparse
import numpy as np
from collections import defaultdict
from collections import Counter
from itertools import combinations

parser = argparse.ArgumentParser(
                    prog='frac_iden_blk')
parser.add_argument('--input', type=str, default="output")
parser.add_argument('--blk_size', type=int, default=500)

args = parser.parse_args()

# open tree sequence
mts = tskit.load(args.input)

G = mts.genotype_matrix()
num_sites, num_samples = G.shape

sites = mts.tables.sites.position.astype(int)
block_idx = sites // args.blk_size
num_blocks = int(mts.sequence_length // args.blk_size)

blocks = [np.flatnonzero(block_idx == b) for b in range(num_blocks)]

num_iden_blk = Counter()
total_blocks = len(blocks)

for block in blocks:
    sub = G[block, :]  
    genotype_bytestrings = [row.tobytes() for row in sub.T]
    iden_groups = defaultdict(list)

    for samp_idx, genotype in enumerate(genotype_bytestrings):
        iden_groups[genotype].append(samp_idx)

    for iden_group in iden_groups.values():
        if len(iden_groups) > 1:
            for i, j in combinations(iden_group, 2):
                num_iden_blk[i, j] += 1

frac_iden_blk = [num_iden_blk.get((i, j), 0) / total_blocks for i,j in combinations(range(num_samples), 2)]

with open(args.input + "_frac_iden_blk", 'wb') as file:
    pickle.dump(frac_iden_blk, file)

print("---fraction of identical blocks computed---")
