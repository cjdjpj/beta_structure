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

# only so no type error, doesn't actually matter - everytime clonal_tmrca is None, frac_clonal is 0
clonal_tmrca = np.array([0 if x is None else x for x in clonal_tmrca])
clonal_pi = clonal_tmrca * params["mu"] * 2

recombinant_tmrca = (tmrca - np.multiply(frac_clonal, clonal_tmrca))/(1-frac_clonal)
recombinant_pi = recombinant_tmrca * params["mu"] * 2

g = sns.jointplot(
    x=clonal_pi,
    y=recombinant_pi,
    kind="scatter",
    height=9, 
    space=0,
    alpha=0,
    marginal_kws={"bins": 160}
)

sns.scatterplot(
    x=clonal_pi,
    y=recombinant_pi,
    hue=frac_clonal,
    palette="viridis",
    ax=g.ax_joint,
    legend=True
)

g.set_axis_labels("Clonal distance", "Recombinant distance", fontsize=12)
g.figure.suptitle("Clonal distance vs recombinant distance (" + run_index + ")", y=1.02)

g.ax_joint.legend(title="Clonal fraction")

if save_fig:
    g.figure.savefig("../figures/" + run_index + "z.png", dpi=300, bbox_inches="tight")
else:
    plt.subplots_adjust(bottom=0.1, left=0.1)
    plt.show()

