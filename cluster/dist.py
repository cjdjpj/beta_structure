import json
import tskit
import pickle
import argparse
from itertools import combinations

parser = argparse.ArgumentParser(
                    prog='dist')
parser.add_argument('--input', type=str, default="output")

args = parser.parse_args()

# open files
mts = tskit.load(args.input)
with open(args.input + ".json", "r") as file:
    params = json.load(file)

# choose pairs
pairs = list(combinations(range(params["nsample"]), 2))

# compute
dist = mts.diversity(pairs, mode='site')

# dump to file
with open(args.input + "_dist", 'wb') as file:
    pickle.dump(dist, file)

print("---distance matrix computed---")
