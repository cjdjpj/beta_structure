import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

header = ["hopkins", "CV"]
df0_003 = pd.read_csv("runs/r_m_tests/r0.003", header = None, names=header)
df0_03 = pd.read_csv("runs/r_m_tests/r0.03", header = None, names=header)
df0_3= pd.read_csv("runs/r_m_tests/r0.3", header = None, names=header)
df3= pd.read_csv("runs/r_m_tests/r3", header = None, names=header)

dfs = [df0_003, df0_03, df0_3, df3]

hopkins = pd.DataFrame({
    "r=0.003": df0_003["hopkins"],
    "r=0.03": df0_03["hopkins"],
    "r=0.3": df0_3["hopkins"],
    "r=3": df3["hopkins"],
})

cv = pd.DataFrame({
    "r=0.003": df0_003["CV"],
    "r=0.03": df0_03["CV"],
    "r=0.3": df0_3["CV"],
    "r=3": df3["CV"],
})

print("---Mean hopkins statistic---")
print(hopkins["r=0.003"].mean())
print(hopkins["r=0.03"].mean())
print(hopkins["r=0.3"].mean())
print(hopkins["r=3"].mean())

print("---Mean CVs---")
print(cv["r=0.003"].mean())
print(cv["r=0.03"].mean())
print(cv["r=0.3"].mean())
print(cv["r=3"].mean())

print("---CV of CVs---")
print(np.std(cv["r=0.003"]) / np.mean(cv["r=0.003"]))
print(np.std(cv["r=0.03"]) / np.mean(cv["r=0.03"]))
print(np.std(cv["r=0.3"]) / np.mean(cv["r=0.3"]))
print(np.std(cv["r=3"]) / np.mean(cv["r=3"]))

plt.figure(figsize = (9,9))
# sns.boxplot(hopkins, showmeans=True, meanprops={'marker':'x','markeredgecolor':'black'})
sns.violinplot(hopkins)
plt.xlabel("Effective recombination rate (r/m)")
plt.ylabel("Hopkins statistic of genetic distance (cluster tendancy)")
sns.pointplot(hopkins, color='black', estimator=np.mean)
plt.show()

plt.figure(figsize = (9,9))
# sns.boxplot(cv, showmeans=True, meanprops={'marker':'x','markeredgecolor':'black'})
sns.violinplot(cv)
plt.xlabel("Effective recombination rate (r/m)")
plt.ylabel("Coefficient of variance of pairwise genetic distance")
sns.pointplot(cv, color='black', estimator=np.mean)
plt.ylim(0, None)
plt.show()
