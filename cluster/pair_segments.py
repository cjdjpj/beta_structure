import random
import pickle
import tskit
import argparse
from itertools import combinations
from dataclasses import dataclass

parser = argparse.ArgumentParser(
                    prog='dist')
parser.add_argument('--input', type=str, default="output")
parser.add_argument('--peak', type=int, default=1)
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

if peak == 1:
    lower_bound = 0.006
    higher_bound = 0.012
elif peak == 2:
    lower_bound = 0.018
    higher_bound = 0.022
elif peak == 3:
    lower_bound = 0.024
    higher_bound = 0.028
else:
    raise RuntimeError("invalid peak, choose one of 1,2,3")

valid_pairs = [i for i, num in enumerate(dist) if lower_bound <= num <= higher_bound]

pairs_of_focus = random.sample(valid_pairs, args.num_pairs)
print("DONE: pair of focus " + str([pairs[p] for p in pairs_of_focus]))

GENE_CONVERSION_FLAG = 1 << 21

gc_nodes = [u.id for u in mts.nodes() if u.flags == GENE_CONVERSION_FLAG]
real_gc = set(gc_nodes[1::2])

class RecombinantType:
    pass

@dataclass
class DoublyRecombined(RecombinantType):
    pass

@dataclass
class SinglyRecombined(RecombinantType):
    sample: int

@dataclass
class NotRecombined(RecombinantType):
    pass

def recomb_seg_status(i, j, real_gc_nodes, tree):
    og_i = i
    og_j = j
    mrca = tree.mrca(i,j)

    # no real gc node on path from i to mrca
    i_recombined = False
    while i != mrca:
        if i in real_gc_nodes:
            i_recombined = True
            break
        i = tree.parent(i)

    # no real gc node on path from j to mrca
    j_recombined = False
    while j != mrca:
        if j in real_gc_nodes:
            j_recombined = True
            break
        j = tree.parent(j)

    if i_recombined and j_recombined:
        return DoublyRecombined()
    elif i_recombined:
        return SinglyRecombined(og_i)
    elif j_recombined:
        return SinglyRecombined(og_j)
    else:
        return NotRecombined()
     
segments_tmrca = []
print("computing tmrcas for pair " + str([pairs[p] for p in pairs_of_focus]))
print("num_trees = ", mts.num_trees) 
k=0
for tree in mts.trees():
    if k%5000 == 0:
        print("at tree = ", k) 
    k+=1
    for i, j in [pairs[p] for p in pairs_of_focus]:
        recomb_status = recomb_seg_status(i, j, real_gc, tree)
        pair_tmrca = tree.tmrca(i,j) 
        if isinstance(recomb_status, NotRecombined):
            segments_tmrca.append(("Clonal", pair_tmrca))
        elif isinstance(recomb_status, SinglyRecombined):
            segments_tmrca.append(("Singly", pair_tmrca, recomb_status.sample))
        else:
            segments_tmrca.append(("Doubly", pair_tmrca))

print("COMPUTATION DONE") 

with open(args.input + "_peak-" + str(peak) + "_pair_segments", 'wb') as file:
    pickle.dump(segments_tmrca, file)
