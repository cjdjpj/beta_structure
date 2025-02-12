import json
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
import numpy as np
from scipy.spatial.distance import squareform
from sklearn.manifold import MDS
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

print(len(distance_list))

print("Average pi:" + str(sum(distance_list)/len(distance_list)))

### PAIRWISE DISTANCE HISTOGRAM
plt.figure(figsize = (9,9))
sns.histplot(distance_list, stat='probability')
plt.xlabel("Pairwise mean number of nucleotide differences (Nei's pi)")
plt.ylabel("Frequency")
if save_fig == True:
    plt.savefig(run_index + "a.png", dpi=300)
else:
    plt.show()

### MDS
mds = MDS(n_components=2, dissimilarity="precomputed", random_state=42)
distance_matrix = squareform(distance_list)
mds_coords = mds.fit_transform(distance_matrix)

plt.figure(figsize = (9,9))
sns.scatterplot(x=mds_coords[:, 0], y=mds_coords[:, 1])
plt.xlabel("MDS 1")
plt.ylabel("MDS 2")
if save_fig == True:
    plt.savefig(run_index + "b.png", dpi=300)
else:
    plt.show()

### PAIRWISE DISTANCE DENDROGRAM
Z = sch.linkage(distance_list, method='average')
mu = params["mu"]
Z[:, 2] = Z[:, 2] / mu

plt.figure(figsize=(9, 9))
sch.dendrogram(Z)
plt.ylabel("Generations")
plt.xlabel("Samples")
plt.xticks([])
if save_fig == True:
    plt.savefig(run_index + "c.png", dpi=300)
else:
    plt.show()

