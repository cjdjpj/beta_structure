import numpy as np
import json
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
np.set_printoptions(legacy='1.25')

###
run_index = "r002"
###

input_path = "runs/" + run_index

with open(input_path + ".json", "r") as file:
    params = json.load(file)
print(json.dumps(params, indent = 4))

with open(input_path + "_frac_clonal", "rb") as file:
    clonal_tmrca = pickle.load(file)

with open(input_path + "_dist", "rb") as file:
    dist = pickle.load(file)

frac_clonal, clonal_tmrca = zip(*clonal_tmrca)
frac_clonal = np.array(frac_clonal)
clonal_tmrca = np.array(clonal_tmrca)

tmrca = dist/(params["mu"]*2)

# only so no type error, doesn't actually matter - everytime clonal_tmrca is None, frac_clonal is 0
clonal_tmrca = [0 if x is None else x for x in clonal_tmrca] 

recombinant_tmrca = (tmrca - np.multiply(frac_clonal, clonal_tmrca))/(1-frac_clonal)
recombinant_pi = recombinant_tmrca * params["mu"] * 2

g = sns.jointplot(
    x=dist,
    y=recombinant_pi,
    kind="scatter",
    height=6, 
    space=0,
    alpha=0,
    marginal_kws={"bins": 160}
)

sns.scatterplot(
    x=dist,
    y=recombinant_pi,
    hue=frac_clonal,
    palette="viridis",
    ax=g.ax_joint,
    legend=True
)

g.set_axis_labels("Pairwise distance", "Recombinant segments distance", fontsize=12)
rho = 2 * params["r"] * params["track_length"] * params["KT_2"]
model_str = "kingman" if params["model"] == "kingman" else "beta ($\\alpha = $" + str(params["alpha"]) + ")" 
g.figure.suptitle("Pairwise distance vs recombinant segments distance (" + model_str + ", $\\rho$=" + str(rho)  + ")")

g.ax_joint.legend(title="Clonal fraction")

# plot y=x
xlim = g.ax_joint.get_xlim()
ylim = g.ax_joint.get_ylim()
lims = [max(min(xlim[0], ylim[0]), -np.inf), min(max(xlim[1], ylim[1]), np.inf)]
g.ax_joint.plot(lims, lims, color='red', linestyle='--', linewidth=1)
g.ax_joint.set_xlim(xlim)
g.ax_joint.set_ylim(ylim)

plt.subplots_adjust(bottom=0.1, left=0.1, top=0.95)
plt.show()

