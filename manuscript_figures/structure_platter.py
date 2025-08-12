import numpy as np
import json
import pickle 
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
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

save_fig = True
run_indices = ["156", "234", "151", "77"]

def load_run(run_index):
    input_path = "runs_structured/" + run_index

    with open(input_path + ".json", "r") as file:
        params = json.load(file)

    with open(input_path + "_dist", "rb") as file:
        dist = pickle.load(file)
        dist = np.array(dist)

    with open(input_path + "_rd", "r") as file:
        r_d = float(file.read())

    with open(input_path + "_frac_clonal", "rb") as file:
        clonal_tmrca = pickle.load(file)

    frac_clonal, clonal_tmrca = map(np.array, zip(*clonal_tmrca))
    clonal_tmrca = np.array([0 if x is None else x for x in clonal_tmrca])

    return dist, frac_clonal, clonal_tmrca, r_d, params

# CREATE FIGURE
fig, axes = plt.subplot_mosaic(
    [
        ["A", "A", "B", "B", "C", "C", "D", "D"],
        ["a", "a", "b", "b", "c", "c", "d", "d"],
        ["a", "a", "b", "b", "c", "c", "d", "d"],
    ],
    figsize = (8, 4),
    sharex = True
)

for label, run_index in zip(["A", "B", "C", "D"], run_indices):
    dist, frac_clonal, clonal_tmrca, r_d, params = load_run(run_index)
    ax = axes[label]

    recomb_status = [
        "Fully recombined" if frac == 0 
        else "Partially\nrecombined" if 0 < frac < 1 
        else "Fully clonal" 
        for frac in frac_clonal
    ]

    sns.histplot(
        x=dist, stat="probability", hue=recomb_status,
        bins=40, multiple="stack", hue_order=["Partially\nrecombined", "Fully recombined", "Fully clonal"],
        ax=ax, legend = False
    )

    ax.text(-0.1, 1.15, rf"$\textbf{{{label}}}$", transform=ax.transAxes, 
            fontweight="bold", va="top", ha="left")
    
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.2f"))
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.1f"))
    ax.set_ylim(0, 0.25)
    ax.set_xlabel("")
    ax.set_ylabel("")
    if label != "A":
        ax.set_yticklabels([])

for label, run_index in zip(["a", "b", "c", "d"], run_indices):
    ax = axes[label]
    dist, frac_clonal, clonal_tmrca, r_d, params = load_run(run_index)

    recomb_status = [
        "Fully recombined" if frac == 0 
        else "Partially\nrecombined" if 0 < frac < 1 
        else "Fully clonal" 
        for frac in frac_clonal
    ]

    tmrca = dist/(params["mu"]*2)
    recombinant_tmrca = (tmrca - np.multiply(frac_clonal, clonal_tmrca))/(1-frac_clonal)
    recombinant_pi = recombinant_tmrca * params["mu"] * 2

    ax.set_ylabel("")
    if label != "a":
        ax.set_yticklabels([])
    ax.set_ylim(0, 0.03)

    sns.scatterplot(x=dist, y=recombinant_pi, hue=recomb_status, ax=ax, 
                    hue_order=["Partially\nrecombined", "Fully recombined", "Fully clonal"],
                    legend= (label == "d"),
                    s=22)

sns.move_legend(axes["d"], "lower left")
axes["a"].set_ylabel("Pairwise differences of recombined regions")
axes["a"].set_yticks([0.00, 0.01, 0.02, 0.03])
fig.text(0.5, 0.07, "Pairwise mean number of nucleotide differences", ha="center")
fig.subplots_adjust(left=0.15, bottom=0.15)

if save_fig:
    plt.savefig("../figures/structure_platter.png", dpi=500, bbox_inches = "tight")
else:
    plt.show()
