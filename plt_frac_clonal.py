import json
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
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
    fraction_clonal = pickle.load(file)

print([f for f in fraction_clonal if f > 0.5])
print("Mean fraction clonal: ", np.mean(fraction_clonal))

### plot
plt.figure(figsize = (9,9))
plt.title("Fraction of clonal descent between pairs of genomes")
plt.xlabel("Proportion of genome")
plt.ylabel("Frequency")
sns.histplot(fraction_clonal, stat='probability')

if save_fig == True:
    plt.savefig(run_index + "e.png", dpi=300)
else:
    plt.show()
