import random, math, re
import json, argparse
import pickle
import numpy as np
import msprime

parser = argparse.ArgumentParser(
                    prog='msprime_simulator')
parser.add_argument('--output', type=str, default="output")
parser.add_argument('--Ne', type=int, default=25000000)
parser.add_argument('--length', type=int, default=5000000)
parser.add_argument('--track_length', type=int, default=5000)
parser.add_argument('--nsample', type=int, default=150)
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
print(r)

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
distance_matrix = np.zeros((nsample, nsample))

distance_matrix[i_indices, j_indices] = distance_list
distance_matrix[j_indices, i_indices] = distance_list

print("distance matrix computed")

### FRACTION OF IDENTICAL BLOCKS
def prop_identical_blk(s1, s2, m):
    matches = 0
    n = len(s1)
    num_blocks = n-m+1

    i=0
    prev_match = False
    while i <= n-m:
        if (s1[i : i+m] == s2[i : i+m]):
            matches += 1
            i += 1
            prev_match = True
        elif prev_match:
            i += m  
        else:
            i += 1

    return matches/num_blocks

liu_and_good_indices = random.sample(range(len(pairs)), 2*int(math.sqrt(len(pairs))))
liu_and_good_pairs = []
genotypes = mts.genotype_matrix()
for k in liu_and_good_indices:
    (i,j) = pairs[k]
    genome_i = ''.join(map(str, genotypes[:, i]))
    genome_j = ''.join(map(str, genotypes[:, j]))

    liu_and_good_pairs.append((prop_identical_blk(genome_i, genome_j, 1000), distance_list[k]))

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
