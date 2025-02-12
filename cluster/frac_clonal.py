import json
import tskit
import pickle
import random
import argparse
from itertools import combinations

parser = argparse.ArgumentParser(
                    prog='frac_clonal')
parser.add_argument('--input', type=str, default="output")
parser.add_argument('--num_pairs', type=int, default=50)

args = parser.parse_args()

# open files
mts = tskit.load(args.input)
with open(args.input + ".json", "r") as file:
    params = json.load(file)

# choose pairs
pairs = list(combinations(range(params["nsample"]), 2))
random.seed(42)
pairs_sample = random.sample(pairs, args.num_pairs)

# compute
clonal_tmrca = []
for pair in pairs_sample:
    subtree = mts.simplify(pair)
    tmrca_values = {}
    for tree in subtree.trees():
        tree_tmrca = tree.tmrca(0, 1)
        tmrca_values[tree_tmrca] = tmrca_values.get(tree_tmrca, 0) + tree.interval.span
    most_common_tmrca = max(tmrca_values, key = tmrca_values.get)
    clonal_interval = tmrca_values[most_common_tmrca]
    clonal_tmrca.append((clonal_interval/params["length"], most_common_tmrca))

# dump to file
with open(args.input + "_frac_clonal", 'wb') as file:
    pickle.dump(clonal_tmrca, file)

print("---fraction clonal computed---")
