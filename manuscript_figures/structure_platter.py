import numpy as np
import json
import pickle 
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from sklearn.manifold import MDS
from scipy.spatial.distance import squareform
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

    with open(input_path + "_rd", "r") as file:
        r_d = float(file.read())

    with open(input_path + "_frac_clonal", "rb") as file:
        clonal_tmrca = pickle.load(file)

    frac_clonal, clonal_tmrca = zip(*clonal_tmrca)
    frac_clonal = np.array(frac_clonal)
    clonal_tmrca = np.array(clonal_tmrca)
    clonal_tmrca = [0 if x is None else x for x in clonal_tmrca]
    clonal_tmrca = np.array(clonal_tmrca)

    recomb_status = [
        "Fully\nrecombined" if frac == 0 
        else "Partially\nrecombined" if 0 < frac < 1 
        else "Fully clonal" 
        for frac in frac_clonal
    ]

    return dist, recomb_status, r_d, params

# CREATE FIGURE
fig, axes = plt.subplot_mosaic(
    [
        ["A", "A", "B", "B", "C", "C", "D", "D"],
        ["A", "A", "B", "B", "C", "C", "D", "D"],
    ],
    figsize = (8, 2),
    sharey = True,
    sharex = True
)

plt.subplots_adjust(wspace=0.3, hspace=0.15)

for ax, run_index, label in zip(axes.values(), run_indices, [r"$\textbf{A}$", r"$\textbf{B}$", r"$\textbf{C}$", r"$\textbf{D}$"]):
    dist, recomb_status, r_d, params = load_run(run_index)

    sns.histplot(
        x=dist, stat="probability", hue=recomb_status,
        bins=40, multiple="stack", hue_order=["Partially\nrecombined", "Fully\nrecombined", "Fully clonal"],
        ax=ax, legend = (ax == axes["D"])
    )

    ax.text(-0.1, 1.1, label, transform=ax.transAxes, 
            fontweight="bold", va="top", ha="left")
    

    if ax == axes["D"]:
        ax.text(0.15, 0.95, f"$\\bar r_d$ = {r_d:.3f}", transform=ax.transAxes,
                verticalalignment="top")
    else:
        ax.text(0.45, 0.95, f"$\\bar r_d$ = {r_d:.3f}", transform=ax.transAxes,
                verticalalignment="top")

    ax.set_xlabel("")

    if ax != axes["A"]:
        ax.set_ylabel("")

    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.2f"))
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.1f"))

sns.move_legend(axes["D"], "center left")

fig.text(0.5, 0.00, "Pairwise mean number of nucleotide differences", ha="center")
fig.subplots_adjust(left=0.15, bottom=0.15)

if save_fig:
    plt.savefig("../figures/structure_platter.png", dpi=500, bbox_inches = "tight")
else:
    plt.show()
