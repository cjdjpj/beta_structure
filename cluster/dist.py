import tskit
import pickle
import argparse
from itertools import combinations

parser = argparse.ArgumentParser(
                    prog='dist')
parser.add_argument('--input', type=str, default="output")

args = parser.parse_args()

# open tree sequence
mts = tskit.load(args.input)

pairs = list(combinations(range(mts.num_samples), 2))

# compute
dist = mts.diversity(pairs, mode='site')

# dump to file
with open(args.input + "_dist", 'wb') as file:
    pickle.dump(dist, file)

print("---distance matrix computed---")
