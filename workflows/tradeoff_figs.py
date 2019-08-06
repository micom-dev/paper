import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import LogFormatterMathtext
import joypy as jp
from scipy.stats import pearsonr

rates = (
    pd.read_csv("data/tradeoff.csv")
    .dropna(subset=["abundance"])
    .rename(columns={"sample": "id", "compartments": "genus"})
)
rates["tradeoff"] = rates.tradeoff.round(2).astype("str")
rates.loc[rates.tradeoff == "nan", "tradeoff"] = "none"

replication = (
    pd.read_csv("data/replication_rates.csv")
    .query("intercept > 1")
    .groupby(["id", "genus"])
    .rate.mean()
    .reset_index()
)
both = pd.merge(rates, replication, on=["id", "genus"])
within_samples = (
    both.groupby(["tradeoff", "id"])
    .apply(
        lambda df: pd.Series(
            pearsonr(df.rate, df.growth_rate) + (df.shape[0],),
            index=["rho", "p", "n"],
        )
    )
    .reset_index()
)
across_samples = (
    both.groupby("tradeoff")
    .apply(
        lambda df: pd.Series(
            pearsonr(df.rate, df.growth_rate) + (df.shape[0],),
            index=["rho", "p", "n"],
        )
    )
    .reset_index()
)

plt.axhline(0, color="k", linestyle="dashed")
g = sns.boxplot(
    "tradeoff",
    "rho",
    data=within_samples[within_samples.n > 4],
    color="w",
    fliersize=0,
)
g = sns.stripplot(
    "tradeoff",
    "rho",
    data=within_samples[within_samples.n > 4],
    color="k",
    jitter=0.2,
    size=3,
    alpha=0.25,
)
plt.xlabel("tradeoff")
plt.ylabel("Pearson rho")
plt.savefig("figures/within_sample_correlations.svg")
plt.close()

plt.axhline(0, color="k", linestyle="dashed")
g = plt.plot(
    "tradeoff", "rho", "ko-", data=across_samples[across_samples.n > 6]
)
plt.xlabel("tradeoff")
plt.ylabel("Pearson rho")
plt.savefig("figures/across_sample_correlations.svg")
plt.close()


rates["log_rates"] = np.log10(rates.growth_rate)
rates.loc[rates.growth_rate <= 0, "log_rates"] = -16

fig, axes = jp.joyplot(
    rates, by="tradeoff", column="log_rates", color="cornflowerblue"
)
lf = LogFormatterMathtext()
xax = axes[-1].xaxis
xax.set_ticklabels(lf(10.0 ** ex) for ex in xax.get_ticklocs())
plt.xlabel("growth rate [1/h]")
plt.ylabel("tradeoff")
plt.savefig("figures/dists.svg")
plt.close()

non_zero = (
    rates.groupby(["id", "tradeoff"])
    .apply(lambda df: (df.growth_rate > 1e-6).sum())
    .reset_index(name="non_zero")
)
n_genus = rates.groupby("id").genus.nunique()
plt.axhline(n_genus.mean(), color="k", linestyle="dashed")
g = sns.boxplot("tradeoff", "non_zero", data=non_zero, color="w", fliersize=0)
g = sns.stripplot(
    "tradeoff",
    "non_zero",
    data=non_zero,
    color="k",
    jitter=0.2,
    size=3,
    alpha=0.25,
)
plt.xlabel("tradeoff")
plt.ylabel("no. genera growing")
plt.savefig("figures/non_zero.svg")
plt.close()

non_zero = (
    rates.groupby(["id", "tradeoff"])
    .apply(lambda df: (df.growth_rate > 1e-6).sum() / df.shape[0])
    .reset_index(name="non_zero")
)
g = sns.boxplot("tradeoff", "non_zero", data=non_zero, color="w", fliersize=0)
g = sns.stripplot(
    "tradeoff",
    "non_zero",
    data=non_zero,
    color="k",
    jitter=0.2,
    size=3,
    alpha=0.25,
)
plt.xlabel("tradeoff")
plt.ylabel("percent genera growing")
plt.savefig("figures/percent_growing.svg")
plt.close()
