import json
import tskit
import pickle
import random
import argparse
from itertools import combinations

parser = argparse.ArgumentParser(
                    prog='dist')
parser.add_argument('--input', type=str, default="output")
parser.add_argument('--num_pairs', type=int, default=0)

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
distance_list = mts.diversity(pairs, mode='site')

# dump to file
with open(args.input + "_dist", 'wb') as file:
    pickle.dump(distance_list, file)

print("---distance matrix computed---")
