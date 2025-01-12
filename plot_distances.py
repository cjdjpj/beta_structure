import json
import tskit
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.manifold import MDS
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform

###
save_fig = False
run_index = "r001"
###

input_path = "runs/" + run_index

with open(input_path + ".json", "r") as file:
    params = json.load(file)
    print(json.dumps(params, indent = 4))

mts = tskit.load(input_path)
nsample = mts.num_samples
with open(input_path + "_ds", "rb") as file:
    distance_list = pickle.load(file)
    distance_matrix = pickle.load(file)

print("Average pi:" + str(sum(distance_list)/len(distance_list)))

### PAIRWISE DISTANCE HISTOGRAM
plt.figure(figsize = (9,9))
sns.histplot(distance_list);
plt.xlabel("Pairwise mean number of nucleotide differences (Nei's pi)")
plt.ylabel("Frequency")
if save_fig == True:
    plt.savefig(run_index + "a.png", dpi=300)
else:
    plt.show()

### MDS
mds = MDS(n_components=2, dissimilarity="precomputed", random_state=42)
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
condensed_distance_matrix = squareform(distance_matrix)
Z = sch.linkage(condensed_distance_matrix, method='ward')

plt.figure(figsize=(9, 9))
sch.dendrogram(Z)
plt.ylabel("Continuous dissimilarity measure using Ward's linkage")
plt.xlabel("Samples")
plt.xticks([])
if save_fig == True:
    plt.savefig(run_index + "d.png", dpi=300)
else:
    plt.show()

