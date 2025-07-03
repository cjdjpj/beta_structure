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

titles = [r"$\rho/\mu = 0.0$", r"$\rho/\mu = 0.01$", r"$\rho/\mu = 0.1$"]

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

for ax, cdf, title, label in zip(axes.values(), cdfs, titles, [r"$\textbf{A}$", r"$\textbf{B}$", r"$\textbf{C}$"]):

    sns.histplot(
        data=cdf, x=cdf.columns[0], bins=40, hue="model", stat="probability",
        multiple="layer", ax=ax, legend = (ax == axes["C"])
    )

    ax.text(-0.1, 1.1, label, transform=ax.transAxes, 
            fontweight="bold", va="top", ha="left")

    ax.set_xlim(-0.01, 0.31)
    ax.set_title(title)
    ax.set_xlabel("")
    if ax != axes["A"]:
        ax.set_ylabel("")


fig.text(0.5, 0.02, "$\\bar r_d$", ha="center")
fig.subplots_adjust(left=0.15, bottom=0.15)

if save_fig:
    plt.savefig("../figures/r_d_distribution.png", dpi=300, bbox_inches = "tight")
else:
    plt.show()
