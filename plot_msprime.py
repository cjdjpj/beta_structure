import json
import tskit
import pickle
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.manifold import MDS
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform

###
save_fig = False
run_index = "r001"
null_index = "r016"
###

input_path = "runs/" + run_index

with open(input_path + ".json", "r") as file:
    params = json.load(file)
    print(json.dumps(params, indent = 4))

mts = tskit.load(input_path)
print(mts)
nsample = mts.num_samples
with open(input_path + "_ds", "rb") as file:
    distance_list = pickle.load(file)
    distance_matrix = pickle.load(file)
    liu_and_good_pairs = pickle.load(file)

print("Average pi:" + str(sum(distance_list)/len(distance_list)))
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

### PAIRWISE DISTANCE HISTOGRAM
plt.figure(figsize = (9,9))
sns.histplot(distance_list, bins = np.linspace(0, 0.04, 100));
plt.xlabel("Pairwise mean number of nucleotide differences (Nei's pi)")
plt.ylabel("Frequency")
if save_fig == True:
    plt.savefig(run_index + "a.png", dpi=300)
else:
    plt.show()

### LIU & GOOD PLOT
# null_path = "runs/" + null_index
# with open(null_path + "_ds", "rb") as file:
#     null_distance_list = pickle.load(file)
#     null_distance_matrix = pickle.load(file)
#     null_liu_and_good_pairs = pickle.load(file)
# null_x_vals, null_y_vals = zip(*null_liu_and_good_pairs)
x_vals, y_vals = zip(*liu_and_good_pairs)

plt.figure(figsize = (9,9))
plt.ylim(0, 0.05)
plt.xlim(0, 1)
sns.scatterplot(y=y_vals, x=x_vals)
# sns.scatterplot(y=null_y_vals, x=null_x_vals, color = "grey")
plt.xlabel("Proportion of 5 sequence blocks identical")
plt.ylabel("Pairwise mean number of nucleotide differences (Nei's pi)")
if save_fig == True:
    plt.savefig(run_index + "c.png", dpi=300)
else:
    plt.show()

# print("cv: " + str(np.std(distance_list) / np.mean(distance_list))) #CoV
