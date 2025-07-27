import numpy as np
import json
import pickle 
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform
from skbio.stats.ordination import pcoa
from sklearn.manifold import TSNE

###
save_fig = False
run_index = "r001"
###

input_path = "runs/" + run_index

with open(input_path + ".json", "r") as file:
    params = json.load(file)
print(json.dumps(params, indent = 4))

with open(input_path + "_dist", "rb") as file:
    dist = pickle.load(file)
with open(input_path + "_rd", "r") as file:
    r_d = float(file.read())

avg_dist = np.mean(dist)
print("Average pi:", avg_dist)

### PAIRWISE DISTANCE HISTOGRAM
plt.figure(figsize = (6,6))
sns.histplot(dist, stat='probability', bins=100)
plt.axvline(x = avg_dist, color = 'red', alpha = 0.3, label = "Average $\\pi$")
plt.xlabel("Pairwise mean number of nucleotide differences (Nei's pi)")
plt.ylabel("Frequency")
plt.text(0.05, 0.8, f"$\\bar r_d$ = {r_d:.3f}", transform=plt.gca().transAxes,
         fontsize=9, verticalalignment='top')
rho = params["r_m"] * params["track_length"] * params["pi"]
model_str = "kingman" if params["model"] == "kingman" else "beta ($\\alpha = $" + str(params["alpha"]) + ")" 
plt.title("Pairwise diversity histogram (" + model_str + ", $\\rho$=" + str(rho)  + ")")
plt.legend()
if save_fig:
    plt.savefig("../figures/runs/" + run_index + "_dist.png", dpi=300)
else:
    plt.show()

# ### ARITY
# import tskit
# mts = tskit.load(input_path)
# max_value = (0, None)
# for tree in mts.trees():
#     arity = [
#         (tree.num_children(node), int(tree.time(node))) for node in tree.nodes() if tree.num_children(node) > 1
#     ]
#     max_candidate = max(arity, key=lambda x: x[0])
#     if max_candidate[0] > max_value[0]:
#         max_value = max_candidate
#         print(max_value)
#
# print(max_value)

# ### PCA
# dist_matrix = squareform(dist)
# pcoa_results = pcoa(dist_matrix)
#
# ### 2D PCA
# pcoa_coords = pcoa_results.samples[['PC1', 'PC2']].values
# variance_explained = pcoa_results.proportion_explained
# pc1_var = variance_explained["PC1"]*100
# pc2_var = variance_explained["PC2"]*100
#
# jitter_strength = 0.0005
# jittered_coords = pcoa_coords + np.random.normal(loc=0, scale=jitter_strength, size=pcoa_coords.shape)
#
# plt.figure(figsize=(6,6))
# sns.scatterplot(x=jittered_coords[:, 0], y=jittered_coords[:, 1])
# plt.xlabel(f"PCA 1 ({pc1_var:.2f}%)")
# plt.ylabel(f"PCA 2 ({pc2_var:.2f}%)")
# plt.title(f"PCoA ({pc1_var+pc2_var:.2f}% variance explained) (" + run_index + ")")
# plt.show()

# ### 3D PCA
# pcoa_coords = pcoa_results.samples[['PC1', 'PC2', 'PC3']].values
# variance_explained = pcoa_results.proportion_explained
# pc1_var = variance_explained["PC1"]*100
# pc2_var = variance_explained["PC2"]*100
# pc3_var = variance_explained["PC3"]*100
#
# fig = plt.figure(figsize=(6, 6))
# ax = fig.add_subplot(111, projection='3d')
#
# ax.scatter(pcoa_coords[:, 0], pcoa_coords[:, 1], pcoa_coords[:, 2],
#            c='dimgray', s=20, alpha=0.3, edgecolors='black', linewidth=0.3)
#
# ax.set_xlabel(f"PC1 ({pc1_var:.2f}%)")
# ax.set_ylabel(f"PC2 ({pc2_var:.2f}%)")
# ax.set_zlabel(f"PC3 ({pc3_var:.2f}%)")
# ax.set_title(f"3D PCoA ({pc1_var+pc2_var+pc3_var:.2f}% variance explained) ({run_index})")
#
# plt.show()

# ### Scree plot
# components = 8
# plt.figure(figsize=(6,6))
# sns.lineplot(
#     x=[f"PC{i+1}" for i in range(components)], 
#     y=variance_explained.iloc[:components], 
# )
# plt.xticks(rotation=90)
# plt.ylim(0,1)
# plt.ylabel("Proportion of variance explained")
# plt.title("Scree plot of PCoA (" + run_index + ")")
# plt.show()
#
# ### PAIRWISE DISTANCE DENDROGRAM
# Z = sch.linkage(dist, method='average')
# mu = params["mu"]
# Z[:, 2] = Z[:, 2] / mu
#
# plt.figure(figsize=(6,6))
# sch.dendrogram(Z)
# plt.ylabel("Generations")
# plt.xlabel("Samples")
# plt.title("Pairwise distance average dendrogram (" + run_index + ")")
# plt.xticks([])
# plt.show()
