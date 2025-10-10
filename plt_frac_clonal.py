import math
import numpy as np
import json
import pickle
import seaborn as sns
import matplotlib.pyplot as plt

###
save_fig = False
run_index = "151"
###

input_path = "runs_structured/" + run_index

with open(input_path + ".json", "r") as file:
    params = json.load(file)
    print(json.dumps(params, indent = 4))

with open(input_path + "_frac_clonal", "rb") as file:
    clonal_tmrca = pickle.load(file)

frac_clonal, clonal_tmrca = zip(*clonal_tmrca)

prop_fully_recombined = sum(x==0 for x in frac_clonal)/math.comb(params["nsample"], 2) * 100
print(f"Fully recombined: {prop_fully_recombined:.2f}%")

g = sns.jointplot(
    x=frac_clonal, 
    y=clonal_tmrca, 
    height=6, 
    space=0,
    xlim=(0,1),
    marginal_kws={"bins": 160}
)

mu = params["mu"]
r = params["r"]
t = params["track_length"]
R = r * t

if R != 0:
    x = np.linspace(1e-3, 1, 100)
    y = np.log(x)/(-2*R)
    g.ax_joint.plot(x, y, color='blue', linestyle = "dashed", alpha = 0.3, label='$e^{-2RT}$')
    g.ax_joint.legend(title="Legend")

## labels
g.set_axis_labels("Proportion of genome", "Generations", fontsize=12)
rho = 2 * params["r"] * params["track_length"] * params["KT_2"]
model_str = "kingman" if params["model"] == "kingman" else "beta ($\\alpha = $" + str(params["alpha"]) + ")" 
g.figure.suptitle("Fraction of clonal interval vs clonal $T_{\\text{mrca}}$ (" + model_str + ", $\\rho$=" + str(rho)  + ")")

if save_fig:
    g.figure.savefig("../figures/runs_full/" + run_index + "_frac_clonal.png", dpi=300, bbox_inches="tight")
else:
    plt.subplots_adjust(bottom=0.1, left=0.1, top=0.95)
    plt.show()
