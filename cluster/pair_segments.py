import random
import pickle
import tskit
import argparse
from itertools import combinations

parser = argparse.ArgumentParser(
                    prog='dist')
parser.add_argument('--input', type=str, default="output")
parser.add_argument('--peak', type=str, default="first_peak")
parser.add_argument('--num_pairs', type=int, default=1)

args = parser.parse_args()

peak = args.peak

with open(args.input + "_dist", "rb") as file:
    dist = pickle.load(file)

print("reading in tree sequence")
mts = tskit.load(args.input)
print("DONE: reading in tree sequence")
n = mts.num_samples

pairs = list(combinations(range(n), 2))

if peak == "first_peak":
    lower_bound = 0.006
    higher_bound = 0.012
elif peak == "second_peak":
    lower_bound = 0.018
    higher_bound = 0.022
elif peak == "third_peak":
    lower_bound = 0.024
    higher_bound = 0.028

valid_pairs = [i for i, num in enumerate(dist) if lower_bound <= num <= higher_bound]

pairs_of_focus = random.sample(valid_pairs, args.num_pairs)
print("DONE: pair of focus " + str([pairs[p] for p in pairs_of_focus]))

GENE_CONVERSION_FLAG = 1 << 21

gc_nodes = [u.id for u in mts.nodes() if u.flags == GENE_CONVERSION_FLAG]
real_gc = set(gc_nodes[1::2])

def are_clonal(i, j, real_gc_nodes, tree):
    mrca = tree.mrca(i,j)

    # no real gc node on path from i to mrca
    while i != mrca:
        if i in real_gc_nodes:
            return False
        i = tree.parent(i)

    # no real gc node on path from j to mrca
    while j != mrca:
        if j in real_gc_nodes:
            return False
        j = tree.parent(j)

    return True
     
segments_tmrca = []
print("computing tmrcas for pair " + str([pairs[p] for p in pairs_of_focus]))
print("num_trees = ", mts.num_trees) 
k=0
for tree in mts.trees():
    if k%5000 == 0:
        print("at tree = ", k) 
    k+=1
    for i, j in [pairs[p] for p in pairs_of_focus]:
        if not are_clonal(i, j, real_gc, tree):
            pair_tmrca = tree.tmrca(i,j)
            segments_tmrca.append(pair_tmrca)

print("COMPUTATION DONE") 

with open(args.input + "_" + peak + "_segments_tmrca", 'wb') as file:
    pickle.dump(segments_tmrca, file)
