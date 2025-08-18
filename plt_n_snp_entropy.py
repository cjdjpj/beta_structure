import math
import pickle
import warnings
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
from collections import Counter
warnings.filterwarnings("ignore")
np.set_printoptions(legacy='1.25')
import scienceplots

plt.style.use("science")
plt.rcParams.update({
    "font.size": 10,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 7,
    "figure.titlesize": 10,
})

###
n_vals = [1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13, 15, 17, 19, 21, 23, 25]
run_index = "runs_structured/149"
nsample = 100
###

with open(run_index + "_entropy", "rb") as file:
    sample_entropy_all_n = pickle.load(file)

plt.figure(figsize = (5,5))
# plt.ylim(0, 10)
plt.xlabel("n in n-SNP")
plt.ylabel("Entropy")
plt.title("n-SNP entropy profile")
for label, values in sample_entropy_all_n.items():
    sns.lineplot(x=n_vals, y=values, label=label)

plt.gca().get_legend().remove()
plt.show()
