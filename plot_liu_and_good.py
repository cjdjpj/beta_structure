import json
import tskit
import pickle
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

###
save_fig = False
run_index = "r007"
###

input_path = "runs/" + run_index
blk_size = 1000 # segment of identity length

with open(input_path + ".json", "r") as file:
    params = json.load(file)
    print(json.dumps(params, indent = 4))

mts = tskit.load(input_path)
nsample = mts.num_samples
with open(input_path + "_ds", "rb") as file:
    distance_list = pickle.load(file)
    distance_matrix = pickle.load(file)
    liu_and_good_pairs = pickle.load(file)

average_divergence = sum(distance_list)/len(distance_list)
print("Average pi:" + str(average_divergence))


### LIU & GOOD PLOT

plt.figure(figsize = (9,9))
plt.ylim(0, 0.05)
plt.xlim(0, 1)

### NULL POINTS
null_index = "r001"
null_path = "runs/" + null_index
with open(null_path + "_ds", "rb") as file:
    null_distance_list = pickle.load(file)
    null_distance_matrix = pickle.load(file)
    null_liu_and_good_pairs = pickle.load(file)
null_x_vals, null_y_vals = zip(*null_liu_and_good_pairs)
sns.scatterplot(y=null_y_vals, x=null_x_vals, color = "grey")

### NULL LINE
n_x = np.linspace(0, 1, 1000)
n_y = -1/blk_size * np.log(n_x)
plt.plot(n_x, n_y, color='grey')

### RECOMBINANT LINE
mu = params["mu"]
r_m = params["r_m"]
l = params["length"]
t = params["track_length"]
R = r_m * mu * t

r_x = np.linspace(0, 1, 1000)
r_y = average_divergence * (1-pow(r_x,R/(R+mu*blk_size)))

plt.plot(r_x, r_y, color='red')

x_vals, y_vals = zip(*liu_and_good_pairs)
sns.scatterplot(y=y_vals, x=x_vals)
plt.xlabel("Proportion of 1kb sequence blocks identical")
plt.ylabel("Pairwise mean number of nucleotide differences (Nei's pi)")
if save_fig == True:
    plt.savefig(run_index + "c.png", dpi=300)
else:
    plt.show()

# print("cv: " + str(np.std(distance_list) / np.mean(distance_list))) #CoV
