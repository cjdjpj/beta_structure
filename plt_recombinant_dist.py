import math
import random
import numpy as np
import json
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
np.set_printoptions(legacy='1.25')

###
save_fig = False
run_index = "r001"
###

input_path = "runs/" + run_index

with open(input_path + ".json", "r") as file:
    params = json.load(file)
print(json.dumps(params, indent = 4))

with open(input_path + "_frac_clonal", "rb") as file:
    clonal_tmrca = pickle.load(file)

with open(input_path + "_dist", "rb") as file:
    dist = pickle.load(file)

frac_clonal, most_common_tmrca = zip(*clonal_tmrca)
frac_clonal = np.array(frac_clonal)
most_common_tmrca = np.array(most_common_tmrca)

random.seed(42)
random_pair_indices = random.sample(range(math.comb(params["nsample"], 2)), len(frac_clonal))

adj_dist = dist[random_pair_indices]
adj_tmrca = adj_dist/(params["mu"]*2)

recombinant_tmrca = (adj_tmrca - np.multiply(frac_clonal, most_common_tmrca))/(1-frac_clonal)

plt.figure(figsize = (9,9))
sns.histplot(recombinant_tmrca, stat='probability', bins=160)
plt.xlabel(f"$T_{{\\text{{mrca}}}}$ of recombined regions")
plt.ylabel("Frequency")
plt.title(f"Recombinant segments $T_{{\\text{{mrca}}}}$ histogram (" + run_index + ")")
if save_fig:
    plt.savefig(run_index + "_recomb_segs", dpi=300)
else:
    plt.show()
