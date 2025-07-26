import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
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

###
file_stem2 = "kingman_0.0"
file_stem3 = "kingman_0.01"
file_stem4 = "kingman_0.1"
file_stem6 = "beta1.1_0.0"
file_stem7 = "beta1.1_0.01"
file_stem8 = "beta1.1_0.1"
###

dir = "runs_mass/"
df2 = pd.read_csv(dir + file_stem2 + ".csv", header=None)
df3 = pd.read_csv(dir + file_stem3 + ".csv", header=None)
df4 = pd.read_csv(dir + file_stem4 + ".csv", header=None)
df6 = pd.read_csv(dir + file_stem6 + ".csv", header=None)
df7 = pd.read_csv(dir + file_stem7 + ".csv", header=None)
df8 = pd.read_csv(dir + file_stem8 + ".csv", header=None)

df2["model"] = r"Kingman"
df3["model"] = r"Kingman"
df4["model"] = r"Kingman"
df6["model"] = r"Beta ($\alpha = 1.1$)"
df7["model"] = r"Beta ($\alpha = 1.1$)"
df8["model"] = r"Beta ($\alpha = 1.1$)"

cdf1 = pd.concat([df2, df6], ignore_index=True)
cdf2 = pd.concat([df3, df7], ignore_index=True)
cdf3 = pd.concat([df4, df8], ignore_index=True)
cdfs = [cdf1, cdf2, cdf3]

titles = [r"$\rho = 0$", r"$\rho = 0.75$", r"$\rho = 7.5$"]

# CREATE FIGURE
fig, axes = plt.subplot_mosaic(
    [
        ["A", "A", "B", "B", "C", "C"],
        ["A", "A", "B", "B", "C", "C"],
        ["a", "a", "b", "b", "c", "c"],
    ],
    figsize = (6, 3),
    sharex = True
)

plt.subplots_adjust(wspace=0.3, hspace=0.15)

for ax, cdf, title, label in zip([axes["A"], axes["B"], axes["C"]], cdfs, titles, [r"$\textbf{A}$", r"$\textbf{B}$", r"$\textbf{C}$"]):

    bins = 30
    sns.histplot(
        data=cdf, x=cdf.columns[0], bins=30, hue="model", stat="probability",
        multiple="layer", ax=ax, legend = (ax == axes["C"])
    )

    ax.text(-0.1, 1.1, label, transform=ax.transAxes, 
            fontweight="bold", va="top", ha="left")

    ax.set_xlim(-0.01, 0.61)
    ax.set_ylim(0.00, 0.53)
    ax.set_title(title)
    ax.set_xlabel("")
    ax.set_ylabel("")
    if ax != axes["A"]:
        ax.set_yticklabels([])

for ax, cdf in zip([axes["a"], axes["b"], axes["c"]], cdfs):

    bins = 30
    sns.histplot(
        data=cdf, x=cdf.columns[0], bins=30, hue="model", stat="probability",
        multiple="layer", ax=ax, legend = False
    )

    ax.set_yscale("log")
    ax.set_xlim(-0.01, 0.61)
    ax.set_ylim(3e-4, 0.99)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis='y', which='major', labelsize=6)
    if ax != axes["a"]:
        ax.set_yticklabels([])


fig.text(0.5, 0.07, "$\\bar r_d$", ha="center")
fig.text(0.09, 0.5, "Probability", va="center", rotation="vertical")
fig.subplots_adjust(left=0.15, bottom=0.15)

if save_fig:
    plt.savefig("../figures/r_d_distribution.png", dpi=500, bbox_inches = "tight")
else:
    plt.show()
