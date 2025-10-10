import json
import pickle
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

###
save_fig = False
run_index = "52"
###

input_path = "runs_structured/" + run_index
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
    ylim=(0,0.065),
    marginal_kws={"bins": 160}
)

## NULL POINTS
sns.scatterplot(x=null_frac_iden_blk, y=null_dist, color='grey')

## NULL LINE
n_x = np.linspace(1e-3, 1, 1000)
n_y = -1/blk_size * np.log(n_x)
g.ax_joint.plot(n_x, n_y, color='grey', linestyle='--')

## RECOMBINANT LINE
def expected_dist(f):
    mu     = params["mu"]
    r      = params["r"]
    t      = params["track_length"]

    # per base rate of replacement by recombination
    R = r * (t) * np.exp(-blk_size/t)

    denom = R + mu * blk_size

    # recombination‐driven divergence + clonal (mutational) divergence
    term_recomb = avg_dist * (1 - f**(R/denom))
    term_mut    = f**(R/denom) * (1 - f**(mu/denom))

    return term_recomb + term_mut

r_x = np.linspace(1e-10, 1, 1000)
r_y = expected_dist(r_x)
g.ax_joint.plot(r_x, r_y, color='red', linestyle='--')

## labels
g.set_axis_labels("Proportion of 1kb base blocks identical", 
                  "Pairwise mean number of nucleotide differences", 
                  fontsize=12)
rho = 2 * params["r"] * params["track_length"] * params["KT_2"]
model_str = "kingman" if params["model"] == "kingman" else "beta ($\\alpha = $" + str(params["alpha"]) + ")" 
g.figure.suptitle("Fraction of identical blocks vs distance (" + model_str + ", $\\rho$=" + str(rho)  + ")")

if save_fig:
    g.figure.savefig("../figures/runs_full/" + run_index + "_frac_iden_blk.png", dpi=300, bbox_inches="tight")
else:
    plt.subplots_adjust(bottom=0.1, left=0.1, top=0.95)
    plt.show()
