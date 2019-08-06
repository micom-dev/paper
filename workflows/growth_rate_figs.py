import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatterMathtext
import seaborn as sns
import joypy as jp
from scipy.stats import pearsonr

rates = (
    pd.read_csv("data/tradeoff.csv")
    .dropna(subset=["abundance"])
    .rename(columns={"sample": "id", "compartments": "genus"})
)

rates["log_rates"] = np.log10(rates.growth_rate)
rates.loc[rates.growth_rate <= 0, "log_rates"] = -16

pos = rates.query("growth_rate > 1e-6 and tradeoff == 0.5")
o = pos.groupby("genus").log_rates.mean().sort_values().index
ax = sns.stripplot(
    "genus", "log_rates", data=pos, order=o, alpha=0.25, jitter=True
)
yax = ax.yaxis
lf = LogFormatterMathtext()
yax.set_ticklabels(lf(10.0 ** ex) for ex in yax.get_ticklocs())
ax.xaxis.set_tick_params(rotation=90)
ax.figure.set_figwidth(12)
ax.figure.set_figheight(5)
sns.despine(left=True, bottom=True)
plt.grid(axis="x", color="gainsboro")
plt.ylabel("growth rate [1/h]")
plt.xlabel("")
plt.tight_layout()
plt.savefig("figures/gcs.svg")
plt.close()

fig = plt.figure()
community_growth = (
    rates[rates.tradeoff == 0.5]
    .groupby("id")
    .apply(lambda df: sum(df.abundance * df.growth_rate))
)
sns.distplot(community_growth)
plt.savefig("figures/community_growth.svg")
plt.close()

fig = plt.figure(figsize=(8, 8))
rs = rates[(rates.growth_rate > 1e-6) & (rates.tradeoff == 0.5)]
g = sns.jointplot(
    np.log10(rs.abundance),
    np.log10(rs.growth_rate),
    kind="kde",
    joint_kws={"n_levels": 100},
)
print(
    "log-rates vs. log-abundance rho=%g (p=%g)."
    % pearsonr(rs.abundance, rs.growth_rate)
)
g.ax_joint.set_xlabel("log abundance")
g.ax_joint.set_ylabel("log growth rate")
plt.savefig("figures/rate_vs_abundance.png", dpi=300)
plt.close()
