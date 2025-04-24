import tskit
import pickle
import argparse
from itertools import combinations

parser = argparse.ArgumentParser(
                    prog='frac_clonal')
parser.add_argument('--input', type=str, default="output")

args = parser.parse_args()

# open tree sequence
mts = tskit.load(args.input)

pairs = list(combinations(range(mts.num_samples), 2))

# compute
clonal_tmrca = []
for (i, j) in pairs:
    tmrca_values = {}
    for tree in mts.trees():
        tree_tmrca = tree.tmrca(i, j)
        tmrca_values[tree_tmrca] = tmrca_values.get(tree_tmrca, 0) + tree.interval.span
    most_common_tmrca = max(tmrca_values, key = tmrca_values.get)
    clonal_interval = tmrca_values[most_common_tmrca]
    clonal_tmrca.append((clonal_interval/mts.sequence_length, most_common_tmrca))
    print(clonal_interval/mts.sequence_length)

# dump to file
with open(args.input + "_frac_clonal", 'wb') as file:
    pickle.dump(clonal_tmrca, file)

print("---fraction clonal computed---")
