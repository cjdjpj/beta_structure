import os
import numpy as np
import json
import pickle 
import seaborn as sns
import matplotlib.pyplot as plt
import scienceplots

plt.style.use('science')

###
save_fig = False
run_index = "149"
###

input_path = "runs_structured/" + run_index

with open(input_path + ".json", "r") as file:
    params = json.load(file)
print(json.dumps(params, indent = 4))

with open(input_path + "_dist", "rb") as file:
    dist = pickle.load(file)
if os.path.exists(input_path + "_rd"):
    with open(input_path + "_rd", "r") as file:
        r_d = float(file.read())

avg_dist = np.mean(dist)
print("Average pi:", avg_dist)

with open(input_path + "_frac_clonal", "rb") as file:
    clonal_tmrca = pickle.load(file)

frac_clonal, clonal_tmrca = zip(*clonal_tmrca)
frac_clonal = np.array(frac_clonal)
clonal_tmrca = np.array(clonal_tmrca)
dist = np.array(dist)

tmrca = dist/(params["mu"]*2)

# only so no type error, doesn't actually matter - everytime clonal_tmrca is None, frac_clonal is 0
clonal_tmrca = [0 if x is None else x for x in clonal_tmrca] 

recomb_status = [
    "Fully recombined" if frac == 0 
    else "Partially recombined" if 0 < frac < 1 
    else "Fully clonal" 
    for frac in frac_clonal
]

### PAIRWISE DISTANCE HISTOGRAM
plt.figure(figsize = (4,4))
sns.histplot(x=dist, stat='probability', hue = recomb_status, bins=100, multiple = "stack", hue_order = ["Partially recombined", "Fully recombined", "Fully clonal"])
plt.xlabel("Pairwise mean number of nucleotide differences")
plt.ylabel("Frequency")
if os.path.exists(input_path + "_rd"):
    plt.text(0.05, 0.75, f"$\\bar r_d$ = {r_d:.3f}", transform=plt.gca().transAxes,
             fontsize=9, verticalalignment='top')
rho = 2 * params["r"] * params["track_length"] * params["KT_2"]
model_str = "kingman" if params["model"] == "kingman" else "beta ($\\alpha = $" + str(params["alpha"]) + ")" 
plt.title("Pairwise diversity histogram (" + model_str + ", $\\rho$=" + str(rho)  + ")")
if save_fig:
    plt.savefig("" + run_index + "_distcolored.png", dpi=300)
else:
    plt.show()
