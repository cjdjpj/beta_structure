import json
import pickle
import math
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
from collections import Counter
from scipy.stats import expon
warnings.filterwarnings("ignore")
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
run_index = "r001"
###

input_path = "runs/" + run_index

with open(input_path + ".json", "r") as file:
    params = json.load(file)
with open(input_path + "_2_snp", "rb") as file:
    two_snps = pickle.load(file)

all_sample_occurances = [num for tup in two_snps for num in tup]

counts = list(Counter(all_sample_occurances).values())
sorted_counts = np.sort(counts)
ccdf_y = 1.0 - np.arange(1, len(sorted_counts)+1) / len(sorted_counts)

# exponential fit
loc, scale = expon.fit(sorted_counts)
x_vals = np.linspace(sorted_counts.min(), sorted_counts.max(), 1000)
ccdf_fit = 1 - expon.cdf(x_vals, loc=loc, scale=scale)

plt.figure(figsize=(5, 5))
sns.scatterplot(x=sorted_counts, y=ccdf_y, edgecolor="none")
plt.plot(x_vals, ccdf_fit, color='red', label="Exponential Fit")

plt.xlabel("2-SNP count")
plt.ylabel("$P(X \\geq x)$")
plt.yscale("log")
plt.xscale("log")
plt.title("Reverse Cumulative Distribution")
plt.legend()
plt.grid(True)
plt.show()

### 2-SNP graph with focal sample ###
def two_snp_neighbours(pairs, target):
    counts = Counter()
    for pair in pairs:
        if target in pair:
            for item in pair:
                if item != target:
                    counts[item] += 1

    return {i: counts.get(i, 0) for i in range(params["nsample"]) if i != target}

for focal_sample in range(20):
    freq_dict = two_snp_neighbours(two_snps, focal_sample)

    G = nx.Graph()

    G.add_node(str(focal_sample))

    for i, (node, freq) in enumerate(freq_dict.items()):
        if node != focal_sample:
            G.add_node(node)
            G.add_edge(str(focal_sample), node, weight=freq)

    pos = {str(focal_sample): (0, 0)}
    for idx, node in enumerate(freq_dict.keys()):
        angle = 2 * math.pi * idx / len(freq_dict)
        pos[node] = (math.cos(angle), math.sin(angle))

    edge_widths = [np.log(d['weight']) for _, _, d in G.edges(data=True)]

    plt.figure(figsize=(9, 9))
    nx.draw(
        G, 
        pos, 
        with_labels=True, 
        width=edge_widths,
        node_size=200,
        node_color="red",
        edge_color="grey",
        font_size=10
    )
    plt.show()
