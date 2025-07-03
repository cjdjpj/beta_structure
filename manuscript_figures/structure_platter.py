import json
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
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

fig, ax = plt.subplot_mosaic(
    [
        ["A", "B"],
        ["C", "D"]
    ],
    figsize = (4.8, 4.8),
    sharey = True,
)
plt.subplots_adjust(wspace=0.15, hspace=0.15)

paths = [
    "_deprecated/runs_full/r035",
    "runs_structured/212",
    "_deprecated/runs_full/r030",
    "runs_structured/151",
]

for ax, path, label in zip(ax.values(), paths, [r"$\textbf{A}$", r"$\textbf{B}$", r"$\textbf{C}$", r"$\textbf{D}$"]):
    with open(path + "_dist", "rb") as file:
        dist = pickle.load(file)
    with open(path + ".json", "r") as file:
        params = json.load(file)

    sns.histplot(
        x=dist, stat="probability", bins=50,
        ax=ax, legend = None
    )
    
    ax.text(-0.1, 1.07, label, transform=ax.transAxes, 
            fontweight="bold", va="top", ha="left")

    ax.text(0.1, 0.93, "$\\alpha = " + str(params["alpha"]) + "$", transform=ax.transAxes, 
            fontweight="bold", va="top", ha="left")
    ax.text(0.1, 0.86, "$\\rho/\\mu= " + str(params["r_m"]) + "$", transform=ax.transAxes, 
            fontweight="bold", va="top", ha="left")

    ax.set_ylabel("")

fig.text(0.01, 0.5, "Probability", va="center", rotation="vertical")
fig.text(0.5, 0.04, "Pairwise mean number of nucleotide differences (Nei's $\\pi$)", ha="center")

if save_fig:
    plt.savefig("../figures/structure_platter.png", dpi=300, bbox_inches = "tight")
else:
    plt.show()
