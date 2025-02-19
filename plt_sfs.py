import numpy as np
import json
import tskit
import matplotlib.pyplot as plt
import seaborn as sns

###
save_fig = False
run_index = "r001"
###

input_path = "runs/" + run_index

with open(input_path + ".json", "r") as file:
    params = json.load(file)
    print(json.dumps(params, indent = 4))

mts = tskit.load(input_path)

afs = mts.allele_frequency_spectrum(mode="site", polarised=True)

afs_values = afs[1:]

plt.figure(figsize=(9, 9))
# plt.ylim(0, 0.1)
sns.barplot(afs)

allele_counts = np.arange(1, len(afs))
expected_sfs = 1 / allele_counts
expected_sfs *= np.sum(afs_values) / np.sum(expected_sfs)  # scale expectation to correct magnitude
plt.plot(allele_counts, expected_sfs, color="red", alpha = 0.4, label="Expected SFS (∝ 1/k)")

plt.xlabel("Derived Allele Count")
plt.ylabel("Number of Variants")
plt.title("Unfolded Allele Frequency Spectrum")
plt.xticks(np.arange(0, len(afs), 50))
plt.show()

