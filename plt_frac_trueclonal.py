import numpy as np
import json, pickle
import seaborn as sns
import matplotlib.pyplot as plt

###
save_fig = False
run_index = "r003"
###

input_path = "runs/" + run_index

with open(input_path + ".json", "r") as file:
    params = json.load(file)
    print(json.dumps(params, indent = 4))

with open(input_path + "_frac_trueclonal", "rb") as file:
    trueclonal_tmrca = pickle.load(file)

frac_trueclonal, most_common_tmrca = zip(*trueclonal_tmrca)

g = sns.jointplot(
    x=frac_trueclonal, 
    y=most_common_tmrca, 
    height=9, 
    space=0,
    xlim=(0,1),
    marginal_kws={"bins": 160}
)

mu = params["mu"]
r_m = params["r_m"]
t = params["track_length"]
R = r_m * mu * t

x = np.linspace(0, 1, 100)
y = np.log(x)/(-2*R)

ax = g.ax_joint
ax.plot(x, y, color='blue', linestyle = "dashed", alpha = 0.3, label=f'$e^{{-2RT}}$')

## labels
g.set_axis_labels("Proportion of genome", "Generations", fontsize=12)
g.figure.suptitle(f"Fraction of true clonal interval vs most common $T_{{\\text{{mrca}}}}$ (" + run_index + ")")
ax.legend(title="Legend")

if save_fig:
    g.figure.savefig("../figures/" + run_index + "g.png", dpi=300, bbox_inches="tight")
else:
    plt.subplots_adjust(bottom=0.1, left=0.1)
    plt.show()
