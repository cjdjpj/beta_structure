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
segments = mts.ibd_segments(store_segments=True)

clonal_interval = []
tmrca = []

for pair, segment_list in segments.items():
    ancestor_nodes = {}
    for segment in segment_list:
        ancestor_nodes[segment.node] = ancestor_nodes.get(segment.node, 0) + (segment.right-segment.left)
    clonal_ancestor = max(ancestor_nodes, key=ancestor_nodes.get)
    clonal_interval.append(ancestor_nodes[clonal_ancestor]/mts.sequence_length)
    tmrca.append(mts.node(clonal_ancestor).time)

print(tmrca)
# dump to file
clonal_tmrca = list(zip(clonal_interval, tmrca))
with open(args.input + "_frac_clonal", 'wb') as file:
    pickle.dump(clonal_tmrca, file)

print("---fraction clonal computed---")
