import numpy as np
import json, pickle
import seaborn as sns
import matplotlib.pyplot as plt

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
avg_dist = np.mean(dist)
avg_tmrca = avg_dist/params["mu"]/2
avg_clonal_tmrca = np.mean(most_common_tmrca)
avg_recombinant_tmrca = avg_tmrca*2 - np.mean(frac_clonal) * avg_clonal_tmrca

plt.figure(figsize = (9,9))
sns.histplot(frac_clonal, stat="probability", bins = 160)
plt.ylabel("Frequency")
plt.xlabel("Fraction of genome clonal")
plt.xlim(-0.03, 1.05)
plt.title("msprime fraction clonal (" + run_index + ")")
if save_fig:
    plt.savefig("../figures/" + run_index + "v.png", dpi=300)
else:
    plt.show()
