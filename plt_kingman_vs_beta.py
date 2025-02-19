import matplotlib.pyplot as plt
import seaborn as sns
import scipy
import msprime
from scipy.stats import linregress

###
save_fig = False
###

def n_beta(a, T2):
    return ((T2 * a * scipy.special.beta(2 - a, a)) / ((1 + 1 / (2**(a - 1) * (a - 1)))**a))**(1 / (a - 1))

model_s = "kingman"
if model_s == "beta":
    alpha = 1.1
    model = msprime.BetaCoalescent(alpha=alpha),
    ne = n_beta(alpha, 1)
    alpha_str = r"($\alpha$ = " + str(alpha) + ")"
else:
    model = None
    alpha = ""
    alpha_str = ""
    ne = 1

ts = msprime.sim_ancestry(
    model = model,
    population_size=ne,
    samples=16, 
    ploidy=1, 
)

# svg_output = ts.draw_svg(size=(2000, 1200), y_axis=True, symbol_size=0, node_labels = {})
# with open(model_s + ".svg", "w") as f:
#     f.write(svg_output)

coalescent_times = [
    tree.time(node) for tree in ts.trees() for node in tree.nodes() if tree.num_children(node) > 1
] 

arity = [
    tree.num_children(node) for tree in ts.trees() for node in tree.nodes() if tree.num_children(node) > 1
]

print(len([time for time in coalescent_times if time < 0.002])/len(coalescent_times))
avg = r"Avg $T_2$: " + "{:.3f}".format(sum(coalescent_times)/len(coalescent_times))

plt.figure(figsize = (9,9))

plt.ylim(0, 1)
plt.xlim(0, 1)

plt.text(0.25, 0.8, avg, fontsize=12)

sns.histplot(coalescent_times, stat="probability", bins = 500)
plt.xlabel("Coalescent Time")
plt.ylabel("Frequency")
plt.title("Distribution of " + model_s + alpha_str + r" coalescent $T_2$")
if save_fig == True:
    plt.savefig(model_s + str(alpha) + ".png", dpi=300)
else:
    plt.show()

slope, intercept, r_value, p_value, std_err = linregress(coalescent_times, arity)
print("r^2: ", r_value)
print("p-value: ", p_value)

# Predict Y values
arity_regress = [slope * x + intercept for x in coalescent_times]

plt.figure(figsize = (9,9))
plt.ylim(0, 200)
plt.xlim(0, 1.0)
sns.scatterplot(x = coalescent_times, y = arity)

plt.plot(coalescent_times, arity_regress, color='red')

plt.show()

