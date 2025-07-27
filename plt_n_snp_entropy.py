import math
import pickle
import warnings
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
from collections import Counter
warnings.filterwarnings("ignore")
np.set_printoptions(legacy='1.25')
import scienceplots

plt.style.use("science")
plt.rcParams.update({
    "font.size": 10,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 7,
    "figure.titlesize": 10,
})

###
n_vals = [1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13, 15, 17, 19, 21, 23, 25]
run_index = "is08"
nsample = 100
###

sample_entropy_all_n = defaultdict(list)

def entropy(tuples_list):
    if not tuples_list:
        return 0.0

    canonical = (tuple(sorted(t)) for t in tuples_list)
    freq = Counter(canonical)
    total = len(tuples_list)
    
    entropy = -sum((count / total) * math.log2(count / total) for count in freq.values())
    return entropy

for n in n_vals:
    input_path = "runs_inf_sites/" + run_index + "_" + str(n) + "_snp"

    with open(input_path, "rb") as file:
        two_snps = pickle.load(file)

    snp_neighbours = defaultdict(list)
    for tup in two_snps:
        for s in tup:
            snp_neighbours[s].append(tup)

    for s in range(nsample):
        if s in snp_neighbours:
            sample_entropy_all_n[s].append(entropy(snp_neighbours[s]))
        else:
            sample_entropy_all_n[s].append(0)


plt.figure(figsize = (5,5))
plt.ylim(0, 10)
plt.xlabel("n in n-SNP")
plt.ylabel("Entropy")
plt.title("n-SNP entropy profile")
for label, values in sample_entropy_all_n.items():
    sns.lineplot(x=n_vals, y=values, label=label)

plt.gca().get_legend().remove()
plt.show()
