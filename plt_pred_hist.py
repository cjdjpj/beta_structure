import numpy as np
import json
import pickle
import seaborn as sns
import matplotlib.pyplot as plt

###
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

# only so no type error, doesn't actually matter - everytime clonal_tmrca is None, frac_trueclonal is 0
clonal_tmrca = [0 if x is None else x for x in clonal_tmrca] 

avg_dist = np.mean(dist)
avg_tmrca = avg_dist/(params["mu"]*2)
avg_clonal_tmrca = np.mean(clonal_tmrca)
avg_recombinant_tmrca = avg_tmrca*2 - np.mean(frac_trueclonal) * avg_clonal_tmrca

predicted_dist = np.multiply(frac_trueclonal, clonal_tmrca)+ (1-frac_trueclonal) * avg_recombinant_tmrca

plt.figure(figsize = (9,9))
sns.histplot(predicted_dist, stat='probability', bins=160)
plt.axvline(x = avg_tmrca, color = 'red', alpha = 0.3, label = "Average $T_{\\text{mrca}}$")
plt.xlabel("Predicted $T_{\\text{mrca}}$")
plt.ylabel("Frequency")
plt.title("Predicted $T_{\\text{mrca}}$ histogram (" + run_index + ")")
plt.legend()
plt.show()
