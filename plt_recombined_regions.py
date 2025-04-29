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

with open(input_path + "_frac_trueclonal", "rb") as file:
    clonal_tmrca = pickle.load(file)

with open(input_path + "_dist", "rb") as file:
    dist = pickle.load(file)

frac_trueclonal, clonal_tmrca = zip(*clonal_tmrca)
frac_trueclonal = np.array(frac_trueclonal)
clonal_tmrca = np.array(clonal_tmrca)

tmrca = dist/(params["mu"]*2)

# only so no type error, doesn't actually matter - everytime clonal_tmrca is None, frac_trueclonal is 0
clonal_tmrca = [0 if x is None else x for x in clonal_tmrca] 

recombinant_tmrca = (tmrca - np.multiply(frac_trueclonal, clonal_tmrca))/(1-frac_trueclonal)

recomb_status = [
    "Fully recombined" if frac == 0 
    else "Partially recombined"
    for frac in frac_trueclonal
]

### PAIRWISE DISTANCE HISTOGRAM
plt.figure(figsize = (9,9))
sns.histplot(x=recombinant_tmrca, stat='probability', hue = recomb_status, bins=160, multiple = "stack", hue_order = ["Partially recombined", "Fully recombined"])
plt.xlabel("$T_{\\text{mrca}}$ of recombined regions")
plt.ylabel("Frequency")
plt.title("Recombinant regions $T_{\\text{mrca}}$ histogram (" + run_index + ")")
if save_fig:
    plt.savefig("../figures/" + run_index + "i.png", dpi=300)
else:
    plt.show()
