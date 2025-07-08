import pickle
import argparse
import tskit
import numpy as np

parser = argparse.ArgumentParser(
                    prog='dist')
parser.add_argument('--input', type=str, default="output")
parser.add_argument('--n', type=int, default=2)

args = parser.parse_args()

# open tree sequence
mts = tskit.load(args.input)

gt = mts.genotype_matrix()

def n_snp(arr, n):
    counts = np.bincount(arr)
    alleles = np.flatnonzero(counts == n)
    return [tuple(np.nonzero(arr == a)[0]) for a in alleles]

n_snps = []
for locus in gt:
    result = n_snp(locus, args.n)
    if result is not None:
        n_snps.extend(result)

# dump to file
with open(args.input + "_" + str(args.n) + "_snp", 'wb') as file:
    pickle.dump(n_snps, file)

print(str(args.n) + "-SNP computed")
