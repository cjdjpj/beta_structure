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
average_tmrca = avg_dist/params["mu"]/2

g = sns.jointplot(
    x=frac_clonal, 
    y=most_common_tmrca, 
    height=9, 
    space=0,
    xlim=(0,1),
    marginal_kws={"bins": 160}
)

plt.axhline(y=average_tmrca, color='red')

## labels
g.set_axis_labels("Proportion of genome", "Generations", fontsize=12)
g.figure.suptitle(f"Fraction of clonal interval vs most common $T_{{\\text{{mrca}}}}$ (" + run_index + ")")

if save_fig == True:
    g.figure.savefig(run_index + "e.png", dpi=300, bbox_inches="tight")
else:
    plt.subplots_adjust(bottom=0.1, left=0.1)
    plt.show()
