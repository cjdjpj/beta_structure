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
run_indices = ["is01", "is02", "is04", "is08"]
###

def entropy(tuples_list):
    if not tuples_list:
        return 0.0

    normalized = [frozenset(t) for t in tuples_list]
    freq = Counter(normalized)
    total = sum(freq.values())

    entropy = -sum((count / total) * math.log2(count / total) for count in freq.values())
    return entropy

fig, axes = plt.subplot_mosaic(
    [
        ["A", "A", "B", "B", "C", "C", "D", "D"],
        ["A", "A", "B", "B", "C", "C", "D", "D"],
    ],
    figsize = (9.6, 2.4),
    sharey = True
)
plt.ylim(0, 10)
for run_index, ax, r_m in zip(run_indices, ["A", "B", "C", "D"], ["0.0", "0.01", "0.1", "0.3"]):
    sample_entropy_all_n = defaultdict(list)

    for n in n_vals:
        input_path = "runs_inf_sites/" + run_index + "_" + str(n) + "_snp"

        with open(input_path, "rb") as file:
            two_snps = pickle.load(file)

        with open(input_path, "rb") as file:
            two_snps = pickle.load(file)

        snp_neighbours = defaultdict(list)
        for tup in two_snps:
            for s in tup:
                snp_neighbours[s].append(tup)

        for s in range(100):
            if s in snp_neighbours:
                sample_entropy_all_n[s].append(entropy(snp_neighbours[s]))
            else:
                sample_entropy_all_n[s].append(0)


    for label, values in sample_entropy_all_n.items():
        sns.lineplot(x=n_vals, y=values, label=label, ax = axes[ax], legend=False)
    axes[ax].set_title("$\\rho/\\mu = $" + r_m)
    axes[ax].set_xlim(0, 26)

fig.text(0.1, 0.5, "Entropy of the n-SNP distribution", va="center", rotation="vertical")
fig.subplots_adjust(left=0.15, bottom=0.15)

plt.show()
