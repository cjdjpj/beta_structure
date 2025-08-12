import math
import pickle
import argparse
from collections import defaultdict
from collections import Counter

###
n_vals = [1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13, 15, 17, 19, 21, 23, 25]
nsample = 100
###

parser = argparse.ArgumentParser(
                    prog='n_snp_entropy')
parser.add_argument('--input', type=str, default="output")

args = parser.parse_args()

sample_entropy_all_n = defaultdict(list)

def entropy(tuples_list):
    if not tuples_list:
        return 0.0

    canonical = (tuple(sorted(t)) for t in tuples_list)
    freq = Counter(canonical)
    total = len(tuples_list)
    
    entropy = -sum((count / total) * math.log2(count / total) for count in freq.values())
    return entropy

for n in n_vals:
    with open(args.input + "_" + str(n) + "_snp", "rb") as file:
        n_snp = pickle.load(file)

    snp_neighbours = defaultdict(list)
    for tup in n_snp:
        for s in tup:
            snp_neighbours[s].append(tup)

    for s in range(nsample):
        if s in snp_neighbours:
            sample_entropy_all_n[s].append(entropy(snp_neighbours[s]))
        else:
            sample_entropy_all_n[s].append(0)


# dump to file
with open(args.input + "_entropy", 'wb') as file:
    pickle.dump(sample_entropy_all_n, file)
