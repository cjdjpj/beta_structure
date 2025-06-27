import scipy
import json
import argparse
import msprime
import numpy as np

parser = argparse.ArgumentParser(
                    prog='sim')
parser.add_argument('--output', type=str, default="output")
parser.add_argument('--length', type=int, default=5000000)
parser.add_argument('--track_length', type=int, default=5000)
parser.add_argument('--nsample', type=int, default=100)
parser.add_argument('--mu', type=float, default=0.025)
parser.add_argument('--r_m', type=float, default=0.00)
parser.add_argument('--model', type=str, default="kingman")
parser.add_argument('--alpha', type=float, default=None)
parser.add_argument('--pi', type=float, default=0.03)
parser.add_argument('--seed', type=int, default=None)
parser.add_argument('--mut_seed', type=int, default=None)
parser.add_argument('--store_gc_nodes', action="store_true")

args = parser.parse_args()

l = args.length  # number of genes
t = args.track_length  # tract length
nsample = args.nsample  # the number of genomes sampled
mu = args.mu  # mutation rate
r_m = args.r_m # r/m

r = r_m * mu # per base recomb rate
print("r  =",r)

def T2(a, N):
    """Returns the expected pairwise coalescence time in a Beta coalescent"""
    return np.power(1 + 1 / np.exp2(a - 1) / (a - 1), a) * np.power(N, a - 1) / a / scipy.special.beta(2 - a, a)

def n_beta(a, T2):
    """Returns the N necessary to make a Beta coalescent with the given alpha = a have the specified pairwise coalescence time"""
    """if returns 0, means N is unimportant. if returns inf, means T unattainable with given alpha"""
    return ((T2 * a * scipy.special.beta(2 - a, a)) / ((1 + 1 / (2**(a - 1) * (a - 1)))**a))**(1 / (a - 1))

if args.model == "kingman":
    model = None
    Ne = args.pi/(2*mu)

elif  args.model == "beta":
    model = msprime.BetaCoalescent(alpha=args.alpha)
    Ne = n_beta(args.alpha, args.pi/(2*mu))

else:
    raise ValueError(f"Invalid model argument: {args.model}")

print("Ne =",Ne)

kwargs = dict(
    samples=nsample,
    model=model,
    population_size=Ne,
    ploidy=1,
    sequence_length=l,
    gene_conversion_rate=r,
    gene_conversion_tract_length=t,
    random_seed = args.seed
)

if args.store_gc_nodes:
    kwargs["additional_nodes"] = (msprime.NodeType.GENE_CONVERSION)
    kwargs["coalescing_segments_only"] = False

ts = msprime.sim_ancestry(**kwargs)

print("---ancestry simulation done---")

mts = msprime.sim_mutations(ts, rate=mu, random_seed = args.mut_seed)

print("---mutation simulation done---")

# save tree_sequence
mts.dump(args.output)

# save params to json
params_dict = vars(args)
with open(args.output + ".json", 'w') as metadata_file:
    json.dump(params_dict, metadata_file, indent=4)
