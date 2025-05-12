from itertools import combinations
import tskit
import numpy as np
import json
import pickle 
import seaborn as sns
import matplotlib.pyplot as plt

###
save_fig = False
run_index = "r015"
###

input_path = "runs/" + run_index

with open(input_path + ".json", "r") as file:
    params = json.load(file)
print(json.dumps(params, indent = 4))

with open(input_path + "_dist", "rb") as file:
    dist = pickle.load(file)
avg_dist = np.mean(dist)

with open(input_path + "_frac_clonal", "rb") as file:
    clonal_tmrca = pickle.load(file)

frac_clonal, clonal_tmrca = zip(*clonal_tmrca)

mts = tskit.load(input_path)

# find first fully recombined pair
for (a, b) in list(combinations(range(mts.num_samples), 2)):
    pair_index = int(b - a - 1 + params["nsample"] * a - (a * (a + 1)) // 2)
    if frac_clonal[pair_index] == 0:
        i = a
        j = b
        print(i)
        print(j)
        break

print(mts.num_trees)
coalescence_times = []
for tree in mts.trees():
    coalescence_times.append(tree.tmrca(i,j))

plt.figure(figsize = (9,9))
sns.histplot(coalescence_times, bins=160)
plt.show()
