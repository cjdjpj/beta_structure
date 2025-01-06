import random, math, re
import json, argparse
import pickle
import numpy as np
import msprime

parser = argparse.ArgumentParser(
                    prog='msprime_simulator')
parser.add_argument('--output', type=str, default="output")
parser.add_argument('--Ne', type=int, default=1)
parser.add_argument('--length', type=int, default=1000)
parser.add_argument('--track_length', type=int, default=1)
parser.add_argument('--nsample', type=int, default=500)
parser.add_argument('--mu', type=float, default=0.0000006)
parser.add_argument('--r_m', type=float, default=0.1)
parser.add_argument('--model', type=str, default="kingman")

args = parser.parse_args()

l = args.length  # number of genes
t = args.track_length  # tract length
nsample = args.nsample  # the number of genomes sampled
mu = args.mu  # mutation rate
r_m = args.r_m # r/m

r = r_m * mu * l

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

mts = msprime.sim_mutations(ts, rate=mu)

pairs = [(i, j) for i in range(nsample) for j in range(i)]
distance_list = mts.diversity(pairs, mode='site')
distance_matrix = np.zeros((nsample,nsample))

for pair_index, distance in enumerate(distance_list):
    (i,j) = pairs[pair_index]
    distance_matrix[i, j] = distance
    distance_matrix[j, i] = distance

## FRACTION OF IDENTICAL BLOCKS
def prop_identical_blk(s1, s2, m):
    matches = 0
    if len(s1)!=len(s2):
        raise Exception("Sequence length doesn't match")
    n = len(s1)
    num_blocks = n // m

    for i in range(num_blocks):
        if s1[i * m : (i + 1) * m] == s2[i * m : (i + 1) * m]:
            matches += 1

    return matches / num_blocks

liu_and_good_indices = random.sample(range(len(pairs)), int(math.sqrt(len(pairs))))
liu_and_good_pairs = []
for k in liu_and_good_indices:
    (i,j) = pairs[k]
    genome_i = next(mts.haplotypes(samples=[i]))
    genome_j = next(mts.haplotypes(samples=[j]))
    
    liu_and_good_pairs.append((prop_identical_blk(genome_i, genome_j, 10), distance_list[k]))

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
