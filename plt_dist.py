import tskit
import json
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
import numpy as np
from scipy.spatial.distance import squareform
from skbio.stats.ordination import pcoa
np.set_printoptions(legacy='1.25')

###
save_fig = False
run_index = "r001"
###

input_path = "runs/" + run_index

with open(input_path + ".json", "r") as file:
    params = json.load(file)
    print(json.dumps(params, indent = 4))

with open(input_path + "_dist", "rb") as file:
    distance_list = pickle.load(file)

mts = tskit.load(input_path)

average_divergence = sum(distance_list)/len(distance_list)
print("Average pi:", average_divergence)

### ARITY
tree = next(mts.trees())
arity = [
    (tree.num_children(node), int(tree.time(node))) for node in tree.nodes() if tree.num_children(node) > 1
]

arity.sort(reverse=True)

print(arity[:(len(arity)//20)])

### PAIRWISE DISTANCE HISTOGRAM
plt.figure(figsize = (9,9))
sns.histplot(distance_list, stat='probability', bins=160)
plt.axvline(x = average_divergence, color = 'red')
plt.xlabel("Pairwise mean number of nucleotide differences (Nei's pi)")
plt.ylabel("Frequency")
plt.title("msprime pairwise diversity histogram (" + run_index + ")")
if save_fig == True:
    plt.savefig(run_index + "a.png", dpi=300)
else:
    plt.show()

### PCA
distance_matrix = squareform(distance_list)
pcoa_results = pcoa(distance_matrix)
pcoa_coords = pcoa_results.samples[['PC1', 'PC2']].values
variance_explained = pcoa_results.proportion_explained
pc1_var = variance_explained["PC1"]*100
pc2_var = variance_explained["PC2"]*100

plt.figure(figsize = (9,9))
sns.scatterplot(x=pcoa_coords[:, 0], y=pcoa_coords[:, 1])
plt.xlabel(f"PCA 1 ({pc1_var:.2f}%)")
plt.ylabel(f"PCA 2 ({pc2_var:.2f}%)")
plt.title(f"PCoA ({pc1_var+pc2_var:.2f}% variance explained) (" + run_index + ")")
if save_fig == True:
    plt.savefig(run_index + "b.png", dpi=300)
else:
    plt.show()

### PCA components
components = 8
plt.figure(figsize = (9,9))
sns.barplot(
    x=[f"PC{i+1}" for i in range(components)], 
    y=variance_explained.iloc[:components], 
)
plt.xticks(rotation=90)
plt.ylim(0,1)
plt.ylabel("Proportion of variance explained")
plt.title("Proportion of variance explained by PCoA (" + run_index + ")")
plt.show()

### PAIRWISE DISTANCE DENDROGRAM
Z = sch.linkage(distance_list, method='average')
mu = params["mu"]
Z[:, 2] = Z[:, 2] / mu

plt.figure(figsize=(9, 9))
sch.dendrogram(Z)
plt.ylabel("Generations")
plt.xlabel("Samples")
plt.title("Pairwise distance average dendrogram (" + run_index + ")")
plt.xticks([])
if save_fig == True:
    plt.savefig(run_index + "c.png", dpi=300)
else:
    plt.show()

