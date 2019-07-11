import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatterMathtext
import seaborn as sns
import joypy as jp
import pickle

community_mass = 200

with open("../results/tradeoff.pickle", "rb") as pick:
    solutions = pickle.load(pick)


def get_members(sol):
    return sol.apply(lambda s: s.members.drop("medium")).iloc[0]


rates = (solutions.groupby(["tradeoff", "sample"]).solution.
         apply(get_members).reset_index())
rates.growth_rate /= community_mass
rates["log_rates"] = np.log10(rates.growth_rate)
rates.loc[rates.growth_rate <= 0, "log_rates"] = -18

fig, axes = jp.joyplot(rates,
                       by="tradeoff", column="log_rates",
                       color="cornflowerblue")
lf = LogFormatterMathtext()
xax = axes[-1].xaxis
xax.set_ticklabels(lf(10.0**ex) for ex in xax.get_ticklocs())
plt.xlabel("growth rate [1/h]")
plt.ylabel("tradeoff")
plt.savefig("dists.svg")
plt.close()

non_zero = rates.groupby(["sample", "tradeoff"]). \
           apply(lambda df: (df.growth_rate > 1e-6).sum()).reset_index(name="non_zero")
n_genus = rates.groupby("sample").compartments.nunique()
plt.axhline(n_genus.mean(), color="k", linestyle="dashed")
g = sns.boxplot("tradeoff", "non_zero", data=non_zero, color="w", fliersize=0)
g = sns.stripplot("tradeoff", "non_zero", data=non_zero, color="k",
                  jitter=0.3, size=3, alpha=0.25)
plt.xlabel("tradeoff")
plt.ylabel("no. genera growing")
plt.savefig("non_zero.svg")
plt.close()

non_zero = rates.groupby(["sample", "tradeoff"]). \
           apply(lambda df: (df.growth_rate > 1e-6).sum() / df.shape[0]). \
           reset_index(name="non_zero")
g = sns.boxplot("tradeoff", "non_zero", data=non_zero, color="w", fliersize=0)
g = sns.stripplot("tradeoff", "non_zero", data=non_zero, color="k",
                  jitter=0.3, size=3, alpha=0.25)
plt.xlabel("tradeoff")
plt.ylabel("percent genera growing")
plt.savefig("percent_growing.svg")
plt.close()

pos = rates.query("growth_rate > 1e-6 and tradeoff == 0.5")
o = pos.groupby("compartments").log_rates.mean().sort_values().index
ax = sns.stripplot("compartments", "log_rates", data=pos,
                   order=o, alpha=0.25, jitter=True)
yax = ax.yaxis
yax.set_ticklabels(lf(10.0**ex) for ex in yax.get_ticklocs())
ax.xaxis.set_tick_params(rotation=90)
ax.figure.set_figwidth(12)
ax.figure.set_figheight(5)
sns.despine(left=True, bottom=True)
plt.grid(axis="x", color="gainsboro")
plt.ylabel("growth rate [1/h]")
plt.xlabel("")
plt.tight_layout()
plt.savefig("gcs.svg")
plt.close()

fig = plt.figure()
community_growth = rates[rates.tradeoff == 0.5].groupby("sample"). \
                   apply(lambda df: sum(df.abundance * df.growth_rate))
sns.distplot(community_growth)
plt.savefig("community_growth.svg")
plt.close()

fig = plt.figure(figsize=(8, 8))
rs = rates[(rates.growth_rate > 1e-6) & (rates.tradeoff != "none")]
g = sns.jointplot(np.log10(rs.abundance), rs.log_rates, kind="kde",
                  joint_kws={"n_levels": 100}, stat_func=None)
g.ax_joint.set_xlabel("log abundance")
g.ax_joint.set_ylabel("log growth rate")
plt.savefig("rate_vs_abundance.png", dpi=300)
plt.close()
