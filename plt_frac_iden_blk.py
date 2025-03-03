import math
import random
import json
import pickle
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

###
save_fig = False
run_index = "r001"
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
print("Average pi:", average_divergence)

random.seed(42)
random_pair_indices = random.sample(range(math.comb(params["nsample"], 2)), len(frac_iden_blk))

### NULL POINTS
null_index = "r001"
null_path = "runs/" + null_index
with open(null_path + "_frac_iden_blk", "rb") as file:
    null_frac_iden_blk = pickle.load(file)
with open(null_path + "_dist", "rb") as file:
    null_distance_list = pickle.load(file)

adjusted_distance_list = distance_list[random_pair_indices]
adjusted_null_distance_list = null_distance_list[random_pair_indices]

## JOINT PLOT
g = sns.jointplot(
    x=frac_iden_blk, 
    y=adjusted_distance_list, 
    height=9, 
    space=0,
    xlim=(0,1),
    ylim=(0,0.05),
    marginal_kws={"bins": 160}
)

## NULL POINTS
sns.scatterplot(x=null_frac_iden_blk, y=adjusted_null_distance_list, color='grey')

## NULL LINE
n_x = np.linspace(1e-10, 1, 1000)
n_y = -1/blk_size * np.log(n_x)
g.ax_joint.plot(n_x, n_y, color='grey', linestyle='--')

## RECOMBINANT LINE
def expected_divergence(frac_iden):
    mu = params["mu"]
    r_m = params["r_m"]
    t = params["track_length"]
    R = r_m * mu * t
    return average_divergence * (1 - pow(frac_iden, R / (R + mu * blk_size)))

r_x = np.linspace(1e-10, 1, 1000)
r_y = expected_divergence(r_x)
g.ax_joint.plot(r_x, r_y, color='red', linestyle='--')

## labels
g.set_axis_labels("Proportion of 1kb sequence blocks identical", 
                  "Pairwise mean number of nucleotide differences (Nei's pi)", 
                  fontsize=12)
g.figure.suptitle(f"Fraction of identical blocks vs divergence ({run_index})")

if save_fig:
    g.figure.savefig(run_index + "d.png", dpi=300, bbox_inches="tight")
else:
    plt.subplots_adjust(bottom=0.1, left=0.1)
    plt.show()
