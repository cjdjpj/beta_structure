import tskit
import argparse
import numpy as np
from itertools import combinations

parser = argparse.ArgumentParser(
                    prog='rd')
parser.add_argument('--input', type=str, default="output")

args = parser.parse_args()

# open tree sequence
mts = tskit.load(args.input)

# genotype matrix 
gt = mts.genotype_matrix()
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

with open(args.input + "_rd", 'w') as file:
    file.write(str(r_d))
