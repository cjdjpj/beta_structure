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

frac_clonal, clonal_tmrca = zip(*clonal_tmrca)
frac_clonal = np.array(frac_clonal)
clonal_tmrca = np.array(clonal_tmrca)

tmrca = dist/(params["mu"]*2)

# only so types work, doesn't actually matter - everytime clonal_tmrca is None, frac_clonal is 0
clonal_tmrca = [0 if x is None else x for x in clonal_tmrca] 

recombinant_tmrca = (tmrca - np.multiply(frac_clonal, clonal_tmrca))/(1-frac_clonal)

plt.figure(figsize = (9,9))
sns.histplot(recombinant_tmrca, stat='probability', bins=160)
plt.xlabel("$T_{\\text{mrca}}$ of recombined regions")
plt.ylabel("Frequency")
plt.title("Recombinant regions $T_{\\text{mrca}}$ histogram (" + run_index + ")")
if save_fig:
    plt.savefig(run_index + "_recomb_regs", dpi=300)
else:
    plt.show()
