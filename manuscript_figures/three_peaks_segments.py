import numpy as np
import pickle 
import matplotlib.pyplot as plt
import seaborn as sns
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

paths = ["151_first_peak_segments_tmrca",
         "151_second_peak_segments_tmrca",
         "151_third_peak_segments_tmrca"]
paths = ["runs_structured/" + p for p in paths]

fig, axes = plt.subplot_mosaic(
    [
        ["A", "A", "B", "B", "C", "C"],
        ["A", "A", "B", "B", "C", "C"],
    ],
    figsize = (6, 2),
    sharey = True,
)

plt.subplots_adjust(wspace=0.55)

recombined_status = ["20\\% recombined", "fully recombined", "fully recombined"]

bin_edges = np.linspace(0, 0.1, 50)
for label, status, input_path in zip(["A", "B", "C"], recombined_status, paths):
    with open(input_path, "rb") as file:
        tmrcas = pickle.load(file)
    run_name = input_path.removeprefix("runs_structured/")

    recomb_type = []
    divergences = []
    for t in tmrcas:
        if t[0] == "Doubly":
            divergences.append(t[1]*2*0.025)
            recomb_type.append(t[0])
        else:
            divergences.append(t[1]*2*0.025)
            recomb_type.append(t[0]+str(t[2]))

    from collections import Counter
    unique_recomb_types = list(Counter(recomb_type).keys())
    sorted_unique_recomb_types = sorted(unique_recomb_types, key=lambda x: (not x.startswith("Doubly"), x))
    sorted_unique_recomb_types.reverse()
    print(sorted_unique_recomb_types)

    color_dict = {
        sorted_unique_recomb_types[0]: sns.color_palette()[2],
        sorted_unique_recomb_types[1]: sns.color_palette()[1],
        sorted_unique_recomb_types[2]: sns.color_palette()[0]
    }

    axes[label].text(0.2, 0.8, status, transform=axes[label].transAxes,
        ha="left", va="top", fontsize=10)
    axes[label].set_ylim(0, 0.75)
    axes[label].set_xlim(-0.004, 0.1)

    sns.histplot(x=divergences, hue=recomb_type, hue_order=sorted_unique_recomb_types, bins=bin_edges, multiple="stack", stat="probability", palette = color_dict, ax = axes[label])
    
    peak = run_name.split("_")[1] + " " + run_name.split("_")[2]
    axes[label].set_title(peak)
    axes[label].set_xticks([0.00, 0.02, 0.04, 0.06, 0.08, 0.1])

# plt.show()
fig.text(0.5, -0.02, "$d$ of recombined trees", ha="center")
plt.savefig("../figures/151_segments_tmrcas.png", dpi=500)
