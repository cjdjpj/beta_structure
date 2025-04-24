import pickle
import seaborn as sns
import matplotlib.pyplot as plt

###
save_fig = False
run_index = "r001"
###

input_path = "runs/" + run_index

with open(input_path + "_frac_clonal", "rb") as file:
    clonal_tmrca = pickle.load(file)

frac_clonal, most_common_tmrca = zip(*clonal_tmrca)

plt.figure(figsize = (9,9))
sns.histplot(frac_clonal, stat="probability", bins = 160)
plt.ylabel("Frequency")
plt.xlabel("Fraction of genome clonal")
plt.xlim(-0.03, 1.05)
plt.title("msprime fraction clonal (" + run_index + ")")
if save_fig:
    plt.savefig("../figures/" + run_index + "v.png", dpi=300)
else:
    plt.show()
