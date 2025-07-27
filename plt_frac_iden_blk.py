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
    dist = pickle.load(file)

avg_dist = np.mean(dist)
print("Average pi:", avg_dist)

### NULL POINTS
null_index = "r001"
null_path = "runs/" + null_index
with open(null_path + "_frac_iden_blk", "rb") as file:
    null_frac_iden_blk = pickle.load(file)
with open(null_path + "_dist", "rb") as file:
    null_dist = pickle.load(file)

## JOINT PLOT
g = sns.jointplot(
    x=frac_iden_blk, 
    y=dist, 
    height=6, 
    space=0,
    xlim=(0,1),
    marginal_kws={"bins": 160}
)

## NULL POINTS
sns.scatterplot(x=null_frac_iden_blk, y=null_dist, color='grey')

## NULL LINE
n_x = np.linspace(1e-3, 1, 1000)
n_y = -1/blk_size * np.log(n_x)
g.ax_joint.plot(n_x, n_y, color='grey', linestyle='--')

## RECOMBINANT LINE
def expected_dist(frac_iden):
    mu = params["mu"]
    r_m = params["r_m"]
    t = params["track_length"]
    R = r_m * mu * t
    return avg_dist * (1 - pow(frac_iden, R / (R + mu * blk_size)))

r_x = np.linspace(1e-10, 1, 1000)
r_y = expected_dist(r_x)
g.ax_joint.plot(r_x, r_y, color='red', linestyle='--')

## labels
g.set_axis_labels("Proportion of 1kb base blocks identical", 
                  "Pairwise mean number of nucleotide differences (Nei's pi)", 
                  fontsize=12)
rho = params["r_m"] * params["track_length"] * params["pi"]
model_str = "kingman" if params["model"] == "kingman" else "beta ($\\alpha = $" + str(params["alpha"]) + ")" 
g.figure.suptitle("Fraction of identical blocks vs distance (" + model_str + ", $\\rho$=" + str(rho)  + ")")

if save_fig:
    g.figure.savefig("../figures/runs_full/" + run_index + "_frac_iden_blk.png", dpi=300, bbox_inches="tight")
else:
    plt.subplots_adjust(bottom=0.1, left=0.1, top=0.95)
    plt.show()
