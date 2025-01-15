import random, math, re
import json, argparse
import pickle
import numpy as np
import msprime
from scipy.spatial.distance import squareform

parser = argparse.ArgumentParser(
                    prog='msprime_simulator')
parser.add_argument('--output', type=str, default="output")
parser.add_argument('--Ne', type=int, default=25000000)
parser.add_argument('--length', type=int, default=5000000)
parser.add_argument('--track_length', type=int, default=5000)
parser.add_argument('--nsample', type=int, default=50)
parser.add_argument('--mu', type=float, default=0.0000000006)
parser.add_argument('--r_m', type=float, default=0.00)
parser.add_argument('--model', type=str, default="kingman")

args = parser.parse_args()

l = args.length  # number of genes
t = args.track_length  # tract length
nsample = args.nsample  # the number of genomes sampled
mu = args.mu  # mutation rate
r_m = args.r_m # r/m

r = r_m * mu
print("r="+str(r))

match = re.match(r"^beta([0-9]*\.[0-9]+)$", args.model)
if match:
    a = float(match.group(1))
    model = msprime.BetaCoalescent(alpha=a)
elif  args.model == "kingman":
    model = None
else:
    raise ValueError(f"Invalid model argument: {args.model}")

ts = msprime.sim_ancestry(nsample,
                          model=model,
                          population_size=args.Ne,
                          ploidy=1,
                          sequence_length=l,
                          gene_conversion_rate=r,
                          gene_conversion_tract_length=t,
                          record_provenance = True,
                          )

print("---ancestry simulation done---")

mts = msprime.sim_mutations(ts, rate=mu)

print("---mutation simulation done---")

i_indices, j_indices = np.tril_indices(nsample, -1)
pairs = np.column_stack((i_indices, j_indices))

distance_list = mts.diversity(pairs, mode='site')
distance_matrix = squareform(distance_list)

print("distance matrix computed")

### FRACTION OF IDENTICAL BLOCKS
blk_size = 1000 # blk_size
sites = [int(site.position) for site in mts.sites()]

# tally number of sites per block
block_indices = np.floor_divide(sites, blk_size)
muts_per_blk = np.bincount(block_indices)
muts_per_blk = muts_per_blk[muts_per_blk != 0]

def prop_identical_blk(s1, s2):
    i = 0
    matches = 0
    for muts in muts_per_blk:
        if np.array_equal(s1[i : i + muts], s2[i : i + muts]):
            matches += 1
        i += muts
    return matches / len(muts_per_blk)

liu_and_good_indices = random.sample(range(len(pairs)), 2*int(math.sqrt(len(pairs))))
liu_and_good_pairs = []
genotypes = mts.genotype_matrix()
for k in liu_and_good_indices:
    (i,j) = pairs[k]
    liu_and_good_pairs.append((prop_identical_blk(genotypes[:, i], genotypes[:, j]), distance_list[k]))

print("liu and good computed")

### save tree_sequence
mts.dump(args.output)

### save expensive data structures
with open(args.output + "_ds", 'wb') as file:
    pickle.dump(distance_list, file)
    pickle.dump(distance_matrix, file)
    pickle.dump(liu_and_good_pairs, file)

### save params to json
params_dict = vars(args)
with open(args.output + ".json", 'w') as metadata_file:
    json.dump(params_dict, metadata_file, indent=4)
