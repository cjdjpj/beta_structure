import scipy
import argparse
import msprime
import numpy as np
from itertools import combinations

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

def r_d(mts, pairs):
    ### COMPUTE R_D (Agapow & Burt 2001)
    gt = mts.genotype_matrix()
    gt = gt[:500, :] # r_d converges quickly, not necessary to use all SNPs
    m, n = gt.shape

    pairs = np.array(list(combinations(range(n), 2)))

    # compute hamming distances counts
    a = gt[:, pairs[:, 0]]
    b = gt[:, pairs[:, 1]]
    dist = (a == b).astype(np.uint8)

    # variance per loci
    variances = np.var(dist, axis=1)

    # sum of covariance of distances
    def sum_covs(data):
        n = data.shape[1]
        means = data.mean(axis=1, keepdims=True)
        Z = data - means                   
        S = Z.sum(axis=0)                     
        total_sum_all = np.dot(S, S)         
        return total_sum_all / n

    sqrt_variances = np.sqrt(variances)
    sum_sqrt_variances = np.sum(sqrt_variances) ** 2

    # r_d (Agapow & Burt 2001)
    r_d = sum_covs(dist)/sum_sqrt_variances

    return r_d

def pi(mts, pairs):
    """COMPUTE pi"""
    dist = mts.diversity(pairs, mode='site')
    pi = np.mean(dist)

    return pi
    
ts = msprime.sim_ancestry(nsample,
                          model=model,
                          population_size=Ne,
                          ploidy=1,
                          sequence_length=l,
                          gene_conversion_rate=r,
                          gene_conversion_tract_length=t,
                          )

pairs = np.array(list(combinations(range(args.nsample), 2)))
mts = msprime.sim_mutations(ts, rate=mu)

with open(args.output, "w") as file:
    file.write(str(pi(mts, pairs)) + "," + str(r_d(mts,pairs)) + "\n")
