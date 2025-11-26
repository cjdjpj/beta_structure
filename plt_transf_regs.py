import numpy as np
import json
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
np.set_printoptions(legacy='1.25')

###
save_fig = False
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
dist = np.array(dist)

tmrca = dist/(params["mu"]*2)

# only so no type error, doesn't actually matter - everytime clonal_tmrca is None, frac_clonal is 0
clonal_tmrca = [0 if x is None else x for x in clonal_tmrca] 

recombinant_tmrca = (tmrca - np.multiply(frac_clonal, clonal_tmrca))/(1-frac_clonal)
recombinant_pi = recombinant_tmrca * params["mu"] * 2

print("Average pi of recomb regions: ", np.mean(recombinant_pi))

recomb_status = [
    "Fully recombined" if frac == 0 
    else "Partially recombined" if 0 < frac < 1 
    else "Fully clonal" 
    for frac in frac_clonal
]

### PAIRWISE DISTANCE HISTOGRAM
plt.figure(figsize = (6,6))
sns.histplot(x=recombinant_pi, stat='probability', bins=160, multiple = "stack", hue = recomb_status, hue_order = ["Partially recombined", "Fully recombined", "Fully clonal"])
plt.xlabel("Diversity of recombined regions")
plt.ylabel("Frequency")
rho = 2 * params["r"] * params["tract_length"] * params["KT_2"]
model_str = "kingman" if params["model"] == "kingman" else "beta ($\\alpha = $" + str(params["alpha"]) + ")" 
plt.title("Diversity of recombinant regions (" + model_str + ", $\\rho$=" + str(rho)  + ")")
if save_fig:
    plt.savefig("../figures/runs_full/" + run_index + "_transf_regs.png", dpi=300)
else:
    plt.show()
