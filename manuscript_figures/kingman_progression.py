import numpy as np
import json
import pickle 
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from scipy.spatial.distance import squareform
from skbio.stats.ordination import pcoa
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
run_indices = ["r001", "r002", "r004"]

def load_run(run_index):
    input_path = "runs/" + run_index

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
        "Fully recombined" if frac == 0 
        else "Partially recombined" if 0 < frac < 1 
        else "Fully clonal" 
        for frac in frac_clonal
    ]

    return dist, recomb_status, r_d, params

# CREATE FIGURE
fig, axes = plt.subplot_mosaic(
    [
        ["A", "A", "B", "B", "C", "C"],
        ["A", "A", "B", "B", "C", "C"],
    ],
    figsize = (7.2, 2.4),
    sharey = True
)

plt.subplots_adjust(wspace=0.3, hspace=0.15)

for ax, run_index, label in zip(axes.values(), run_indices, [r"$\textbf{A}$", r"$\textbf{B}$", r"$\textbf{C}$"]):
    dist, recomb_status, r_d, params = load_run(run_index)

    sns.histplot(
        x=dist, stat="probability", hue=recomb_status,
        bins=40, multiple="stack", hue_order=["Partially recombined", "Fully recombined", "Fully clonal"],
        ax=ax, legend = (ax == axes["C"])
    )

    if ax == axes["A"]:
        inset_ax = ax.inset_axes([0.05, 0.45, 0.35, 0.35])

        dist_matrix = squareform(dist)
        pcoa_results = pcoa(dist_matrix)

        ### 2D PCA
        pcoa_coords = pcoa_results.samples[["PC1", "PC2"]].values
        variance_explained = pcoa_results.proportion_explained

        jitter_strength = 0.001
        jittered_coords = pcoa_coords + np.random.normal(loc=0, scale=jitter_strength, size=pcoa_coords.shape)

        sns.scatterplot(x=jittered_coords[:, 0], y=jittered_coords[:, 1],
                        ax = inset_ax, color = sns.color_palette()[2],
                        edgecolor="black",
                        linewidth = 0.2,
                        s = 15)

        inset_ax.set_xticklabels([])
        inset_ax.set_yticklabels([])
        inset_ax.set_xlabel("")
        inset_ax.set_ylabel("")

    ax.text(-0.1, 1.1, label, transform=ax.transAxes, 
            fontweight="bold", va="top", ha="left")
    
    ax.text(0.05, 0.95, f"$\\bar r_d$ = {r_d:.3f}", transform=ax.transAxes,
            verticalalignment="top")
    
    ax.set_title(f"$\\rho/\\mu = {round(params["r_m"], 3)}$")
    ax.set_xlabel("")

    if ax != axes["A"]:
        ax.set_ylabel("")

    formatter = mticker.FormatStrFormatter("%.2f")
    ax.xaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(formatter)

sns.move_legend(axes["C"], "center left")

fig.text(0.5, 0.02, "Pairwise mean number of nucleotide differences (Nei's $\\pi$)", ha="center")
fig.subplots_adjust(left=0.15, bottom=0.15)

if save_fig:
    plt.savefig("../figures/kingman_progression.png", dpi=300, bbox_inches = "tight")
else:
    plt.show()
