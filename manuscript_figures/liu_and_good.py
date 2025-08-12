import json
import pickle
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
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
save_fig = True
kingman_indices = ["r001", "r002", "r003", "r004", "r008"]
beta_index = "151"
blk_size = 1000 # for analytical prediction
###

## RECOMBINANT LINE
def expected_dist(f):
    mu     = params["mu"]
    r_m    = params["r_m"]
    t      = params["track_length"]

    # per base rate of replacement by recombination
    R = r_m * mu * t

    denom = R + mu * blk_size

    term_recomb = avg_dist * (1 - f**(R/denom)) # SNPs introduced by recombination
    term_mut    = f**(R/denom) * (1 - f**(mu/denom)) # SNPs introduced by mutation

    return term_recomb + term_mut

# CREATE FIGURE
fig, axes = plt.subplot_mosaic(
    [
        ["A", "A", "B", "B", "B", "C"],
        ["A", "A", "B", "B", "B", "C"],
        ["A", "A", "B", "B", "B", "C"]
    ],
    figsize = (6, 3),
    sharey = True,
)

plt.subplots_adjust(wspace=0.35)

# SUBPLOT E: KC
axes["A"].set_xlim(0, 1)
axes["A"].text(-0.12, 1.041, r"$\textbf{A}$", transform=axes["A"].transAxes, 
        fontweight="bold", va="top", ha="left")
axes["A"].set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
fig.text(0.24, 0.91, "KC", ha="center")

for kingman_index, color in zip(kingman_indices, [sns.color_palette()[6], sns.color_palette()[1], sns.color_palette()[2], sns.color_palette()[3], sns.color_palette()[4]]):
    input_path = "runs/" + kingman_index
    with open(input_path + ".json", "r") as file:
        params = json.load(file)
    with open(input_path + "_frac_iden_blk", "rb") as file:
        frac_iden_blk = pickle.load(file)
    with open(input_path + "_dist", "rb") as file:
        dist = pickle.load(file)
    avg_dist = np.mean(dist)

    # simulation points
    rho = params["r_m"] * params["track_length"] * params["pi"]
    sns.scatterplot(x=frac_iden_blk, y=dist,
                    ax=axes["A"],
                    s=8, 
                    color = color,
                    label=f"$\\rho$ = {rho:.3g}")

    # expectation lines
    r_x = np.linspace(1e-10, 1, 1000)
    r_y = expected_dist(r_x)
    axes["A"].plot(r_x, r_y, 
                   linestyle='--', 
                   color=color, 
                   path_effects=[pe.Stroke(linewidth=1.5, foreground='white'), pe.Normal()])

# legend adjustments
sns.move_legend(axes["A"], "upper left")
handles, labels = axes["A"].get_legend_handles_labels()
new_handles = [
    Line2D(
        [0], [0],
        marker='o',
        color='w',
        label=label,
        markerfacecolor=handle.get_facecolor()[0],
        markersize=(28**0.5),
        linestyle='None'
    )
    for handle, label in zip(handles, labels)
]
axes["A"].legend(handles=new_handles, labels=labels)

### SUBPLOT A - B: BC
input_path = "runs_structured/" + beta_index

with open(input_path + ".json", "r") as file:
    params = json.load(file)
with open(input_path + "_frac_iden_blk", "rb") as file:
    frac_iden_blk = pickle.load(file)
with open(input_path + "_dist", "rb") as file:
    dist = pickle.load(file)
avg_dist = np.mean(dist)

# BC points
sns.scatterplot(x=frac_iden_blk, y=dist, ax=axes["B"], s=8)

# expectation line
r_x = np.linspace(1e-10, 1, 1000)
r_y = expected_dist(r_x)
axes["B"].plot(r_x, r_y, 
               linestyle='--', 
               path_effects=[pe.Stroke(linewidth=1.5, foreground='white'), pe.Normal()])

# BC marginal histogram
sns.histplot(y = dist, bins=60, ax = axes["C"], stat="probability")

# settings
axes["B"].set_yticks([0.0, 0.01, 0.02, 0.03, 0.04])
axes["B"].set_xlim(0, 1)
axes["B"].set_ylim(0, 0.042)
axes["B"].text(-0.08, 1.04, r"$\textbf{B}$", transform=axes["B"].transAxes, 
        fontweight="bold", va="top", ha="left")
axes["C"].set_xlabel("")

alpha = params["alpha"]
fig.text(0.58, 0.91, f"BC ($\\alpha = {alpha:.3g}, \\rho = {rho:.3g}$)", ha="center")
fig.text(0.4, -0.01, "Fraction of 1 kb blocks identical", ha="center")
fig.text(0.04, 0.5, "Pairwise mean number of nucleotide differences", rotation="vertical", va="center")

if save_fig:
    plt.savefig("../figures/liu_and_good.png", dpi=500, bbox_inches = "tight")
else:
    plt.show()
