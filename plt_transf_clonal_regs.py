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

frac_clonal, clonal_tmrca = zip(*clonal_tmrca)
frac_clonal = np.array(frac_clonal)
clonal_tmrca = np.array(clonal_tmrca)

tmrca = dist/(params["mu"]*2)

# remove nones
no_nones = np.array([i for i, val in enumerate(clonal_tmrca) if val is not None])
clonal_tmrca = np.array([clonal_tmrca[i] for i in no_nones])
tmrca = np.array([tmrca[i] for i in no_nones])
frac_clonal = np.array([frac_clonal[i] for i in no_nones])

clonal_pi = clonal_tmrca * params["mu"] * 2
recombinant_tmrca = (tmrca - np.multiply(frac_clonal, clonal_tmrca))/(1-frac_clonal)
recombinant_pi = recombinant_tmrca * params["mu"] * 2

g = sns.jointplot(
    x=clonal_tmrca,
    y=recombinant_tmrca,
    kind="scatter",
    height=6, 
    space=0,
    alpha=0,
    marginal_kws={"bins": 160}
)

sns.scatterplot(
    x=clonal_tmrca,
    y=recombinant_tmrca,
    hue=frac_clonal,
    palette="viridis",
    ax=g.ax_joint,
    legend=True
)

g.set_axis_labels("Clonal $T_{\\text{mrca}}$", "Recombinant $T_{\\text{mrca}}$", fontsize=12)
rho = 2*params["pi"] * params["r_m"]
model_str = "kingman" if params["model"] == "kingman" else "beta ($\\alpha = $" + str(params["alpha"]) + ")" 
g.figure.suptitle("Clonal $T_{\\text{mrca}}$ vs recombinant regions $T_{\\text{mrca}}$ (" + model_str + ", $\\rho$=" + str(rho)  + ")")

g.ax_joint.legend(title="Clonal fraction")

if save_fig:
    g.figure.savefig("../figures/" + run_index + "z.png", dpi=300, bbox_inches="tight")
else:
    plt.subplots_adjust(bottom=0.1, left=0.1, top=0.95)
    plt.show()

