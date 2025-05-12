import math
import tskit
import pickle
import argparse
import numpy as np
import rustworkx as rx
from itertools import combinations

parser = argparse.ArgumentParser(
                    prog='frac_clonal')
parser.add_argument('--input', type=str, default="output")

args = parser.parse_args()

# open tree sequence
mts = tskit.load(args.input)

n = mts.num_samples
L = mts.sequence_length
GENE_CONVERSION_FLAG = 1 << 21

gc_nodes = [u.id for u in mts.nodes() if u.flags == GENE_CONVERSION_FLAG]
real_gc = set(gc_nodes[1::2])

num_pairs = math.comb(n, 2)
clonal_interval = np.zeros(num_pairs, dtype=float)
tmrca = np.full(num_pairs, None, dtype=object)

# for each marginal tree, find components when you "remove" real_gc nodes
for tree in mts.trees():
    g = rx.PyGraph()

    node_map = {} # get graph idx for node idx

    for u in tree.nodes():
        if u not in real_gc:
            node_index = g.add_node(u)
            node_map[u] = node_index

    for u in node_map:
        p = tree.parent(u)
        if p != tskit.NULL and p in node_map:
            g.add_edge(node_map[u], node_map[p], None)

    components = rx.connected_components(g)

    # for all pairs in component, add tree span as clonal
    for comp in components:
        sample_nodes = [g[i] for i in comp if tree.is_sample(g[i])]
        for i, j in combinations(sorted(sample_nodes), 2):
            pair_index = int(j - i - 1 + n * i - (i * (i + 1)) // 2)
            clonal_interval[pair_index] += tree.span
            tmrca[pair_index] = tree.tmrca(i,j)


clonal_tmrca = list(zip(clonal_interval/L, tmrca))

with open(args.input + "_frac_clonal", 'wb') as file:
    pickle.dump(clonal_tmrca, file)

print("---fraction clonal computed---")
