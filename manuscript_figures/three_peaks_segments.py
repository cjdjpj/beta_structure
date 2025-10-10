import numpy as np
import pickle
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
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
np.set_printoptions(legacy="1.25")

save_fig = True

paths = [
    "runs_structured/151_peak-1_pair_segments",
    "runs_structured/151_peak-2_pair_segments",
    "runs_structured/151_peak-3_pair_segments",
]
peaks = [
    "Fig. 3B left peak",
    "Fig. 3B middle peak",
    "Fig. 3B right peak",
]

fig, axes = plt.subplot_mosaic(
    [["D", "D", "E", "E", "F", "F"],
     ["D", "D", "E", "E", "F", "F"]],
    figsize=(7, 2), sharey=True
)
plt.subplots_adjust(wspace=0.65)

bin_edges = np.linspace(0, 0.1, 30)

for label, peak, input_path in zip(["D", "E", "F"], peaks, paths):
    with open(input_path, "rb") as f:
        tmrcas = pickle.load(f)

    divergences = [t[1]*2*0.025 for t in tmrcas]
    recomb_types = []
    singly_clonal_types = {}
    for t in tmrcas:
        if t[0] in ["Doubly", "Clonal"]:
            recomb_types.append(t[0])
        else:
            if t[2] in singly_clonal_types:
                recomb_types.append(t[0] + " " + singly_clonal_types[t[2]])
            else:
                singly_clonal_types[t[2]] = "A" if len(singly_clonal_types) == 0 else "B"
                recomb_types.append(t[0] + " " + singly_clonal_types[t[2]])

    unique_recomb_types = set(recomb_types)

    order = ["Clonal", "Singly B", "Singly A", "Doubly"]

    unique_types = [r for r in order if r in recomb_types]

    color_dict = {
        unique_types[-3]: sns.color_palette()[2],
        unique_types[-2]: sns.color_palette()[1],
        unique_types[-1]: sns.color_palette()[0],
    }

    if len(unique_types) > 3:
        color_dict[unique_types[0]] = sns.color_palette()[3]

    ax = axes[label]
    sns.histplot(x=divergences, hue=recomb_types, hue_order=unique_types,
                 bins=bin_edges, multiple="stack", stat="probability",
                 palette=color_dict, legend=False, ax=ax)

    if label == "D":
        color_dict = {
            unique_types[-3]: (1.0, 0.68627451, 0.25098039),
            unique_types[-2]: (0.25098039, 0.79215686, 0.45490196),
            unique_types[-1]: (0.28627451, 0.5254902, 0.73333333),
            unique_types[0]:  (1.0, 0.37647059, 0.25098039),
        }

        handles = [mpatches.Patch(color=color_dict[cat], label=cat) for cat in unique_types]

        leg = ax.legend(handles=handles, labels=["Clonal", "Singly recombined", "Singly recombined", "Doubly recombined"], loc="upper right")

        for legend_handle in leg.legend_handles:
            legend_handle.set_edgecolor("black")
            legend_handle.set_linewidth(1.1)

    # panel labels
    ax.text(-0.1, 1.1, rf"$\textbf{{{label}}}$", transform=ax.transAxes, 
            fontweight="bold", va="top", ha="left")

    ax.set_title(peak)
    ax.set_ylim(0, 0.82)
    ax.set_xlim(-0.004, 0.1)
    ax.set_xticks(np.arange(0, 0.11, 0.03))

fig.text(0.5, -0.04, "$d$ of segments for all pairs", ha="center")
if save_fig:
    plt.savefig("../figures/151_segments_tmrcas.png", dpi=500)
else:
    plt.show()
