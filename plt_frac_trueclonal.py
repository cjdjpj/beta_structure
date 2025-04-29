import math
import numpy as np
import json
import pickle
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

with open(input_path + "_frac_trueclonal", "rb") as file:
    trueclonal_tmrca = pickle.load(file)

frac_trueclonal, clonal_tmrca = zip(*trueclonal_tmrca)

prop_fully_recombined = sum(x==0 for x in frac_trueclonal)/math.comb(params["nsample"], 2) * 100
print(f"Fully recombined: {prop_fully_recombined:.2f}%")

g = sns.jointplot(
    x=frac_trueclonal, 
    y=clonal_tmrca, 
    height=9, 
    space=0,
    xlim=(0,1),
    marginal_kws={"bins": 160}
)

mu = params["mu"]
r_m = params["r_m"]
t = params["track_length"]
R = r_m * mu * t

if R != 0:
    x = np.linspace(1e-3, 1, 100)
    y = np.log(x)/(-2*R)
    g.ax_joint.plot(x, y, color='blue', linestyle = "dashed", alpha = 0.3, label='$e^{-2RT}$')
    g.ax_joint.legend(title="Legend")

## labels
g.set_axis_labels("Proportion of genome", "Generations", fontsize=12)
g.figure.suptitle("Fraction of true clonal interval vs clonal $T_{\\text{mrca}}$ (" + run_index + ")")

if save_fig:
    g.figure.savefig("../figures/" + run_index + "h.png", dpi=300, bbox_inches="tight")
else:
    plt.subplots_adjust(bottom=0.1, left=0.1)
    plt.show()

# ### plot frac_trueclonal histogram alone
# plt.figure(figsize = (9,9))
# sns.histplot(frac_trueclonal, stat="probability", bins = 160)
# plt.ylabel("Frequency")
# plt.xlabel("Fraction of genome clonal")
# plt.xlim(-0.03, 1.05)
# plt.title("Fraction inferred clonal (" + run_index + ")")
# plt.show()
