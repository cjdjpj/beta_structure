import random
import json
import pickle
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import combinations

###
save_fig = False
run_index = "r002"
###

input_path = "runs/" + run_index
blk_size = 1000 # for analytical prediction

with open(input_path + ".json", "r") as file:
    params = json.load(file)
    print(json.dumps(params, indent = 4))

with open(input_path + "_frac_iden_blk", "rb") as file:
    frac_iden_blk = pickle.load(file)

with open(input_path + "_dist", "rb") as file:
    distance_list = pickle.load(file)

average_divergence = sum(distance_list)/len(distance_list)
print("Average pi:" + str(average_divergence))


### LIU & GOOD PLOT
plt.figure(figsize = (9,9))
plt.ylim(0, 0.05)
plt.xlim(0, 1)

### NULL POINTS
null_index = "r001"
null_path = "runs/" + null_index
with open(null_path + "_frac_iden_blk", "rb") as file:
    null_frac_iden_blk = pickle.load(file)
null_x_vals, null_y_vals = zip(*null_frac_iden_blk)
sns.scatterplot(y=null_y_vals, x=null_x_vals, color = "grey")

### NULL LINE
n_x = np.linspace(0, 1, 1000)
n_y = -1/blk_size * np.log(n_x)
plt.plot(n_x, n_y, color='grey')

### RECOMBINANT LINE
mu = params["mu"]
r_m = params["r_m"]
l = params["length"]
t = params["track_length"]
R = r_m * mu * t

r_x = np.linspace(0, 1, 1000)
r_y = average_divergence * (1-pow(r_x,R/(R+mu*blk_size)))

plt.plot(r_x, r_y, color='red')

pairs = list(combinations(range(params["nsample"]), 2))
random.seed(42)
random_pair_indices = random.sample(range(len(pairs)), len(frac_iden_blk))
adjusted_distance_list = []
for i in random_pair_indices:
    adjusted_distance_list.append(distance_list[i])

sns.scatterplot(y=adjusted_distance_list, x=frac_iden_blk)
plt.xlabel("Proportion of 1kb sequence blocks identical")
plt.ylabel("Pairwise mean number of nucleotide differences (Nei's pi)")
if save_fig == True:
    plt.savefig(run_index + "d.png", dpi=300)
else:
    plt.show()
