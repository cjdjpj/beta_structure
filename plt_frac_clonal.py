import json
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
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
    distance_list = pickle.load(file)

frac_clonal, most_common_tmrca = zip(*clonal_tmrca)
frac_clonal = np.array(frac_clonal)
most_common_tmrca = np.array(most_common_tmrca)

### clonal_interval histogram
plt.figure(figsize = (9,9))
plt.title("Fraction of clonal descent between pairs of genomes (" + run_index + ")")
plt.xlabel("Proportion of genome")
plt.ylabel("Frequency")
sns.histplot(frac_clonal, stat='probability')

if save_fig == True:
    plt.savefig(run_index + "e.png", dpi=300)
else:
    plt.show()

plt.figure(figsize = (9,9))
plt.xlim(0,1)
plt.title(f"Fraction of clonal interval vs most common $T_{{\\text{{mrca}}}}$ (" + run_index + ")")
plt.xlabel("Proportion of genome")
plt.ylabel("Generations")
sns.scatterplot(x=frac_clonal, y=most_common_tmrca)

if save_fig == True:
    plt.savefig(run_index + "f.png", dpi=300)
else:
    plt.show()
