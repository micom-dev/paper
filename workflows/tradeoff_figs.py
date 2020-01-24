import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from plotnine import *
from mizani.formatters import percent_format
from matplotlib.ticker import LogFormatterMathtext
import joypy as jp
from scipy.stats import spearmanr, pearsonr

theme_set(theme_minimal())

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


def corr(df):
    """Calculates correlations.

    In order to run the analysis with Spearman correlations just substitute
    `spearmanr` for `pearsonr` here.
    """
    if df.shape[0] > 2:
        return pd.Series(
            pearsonr(df.rate, df.growth_rate) + (df.shape[0],),
            index=["rho", "p", "n"],
        )
    else:
        return pd.Series(
            [np.nan, np.nan, df.shape[0]], index=["rho", "p", "n"]
        )


both = pd.merge(rates, replication, on=["id", "genus"])
within_samples = both.groupby(["tradeoff", "id"]).apply(corr).reset_index()
across_samples = both.groupby("tradeoff").apply(corr).reset_index()

pl = (
    ggplot(within_samples[within_samples.n > 6], aes(x="tradeoff", y="rho"))
    + geom_hline(yintercept=0, linetype="dashed")
    + geom_boxplot(outlier_color="none")
    + geom_jitter(width=0.15, height=0, alpha=0.5, stroke=0)
    + labs(x="tradeoff", y="Pearson rho")
)
pl.save("figures/within_sample_correlations.svg", width=5, height=5)
within_samples.to_csv("data/correlation_per_sample.csv")

pl = (
    ggplot(across_samples[across_samples.n > 6], aes(x="tradeoff", y="rho"))
    + geom_hline(yintercept=0, linetype="dashed")
    + geom_point()
    + geom_line(aes(group=1))
    + labs(x="tradeoff", y="Pearson rho")
)
pl.save("figures/across_sample_correlations.svg", width=5, height=5)
across_samples.to_csv("data/correlation_all.csv")

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
    .apply(lambda df: (df.growth_rate > 1e-6).sum() / df.shape[0])
    .reset_index(name="non_zero")
)

pl = (
    ggplot(non_zero, aes(x="tradeoff", y="non_zero"))
    + geom_boxplot(outlier_color="none")
    + geom_jitter(width=0.15, height=0, alpha=0.5, stroke=0)
    + scale_y_continuous(labels=percent_format())
    + labs(x="tradeoff", y="percent taxa growing")
)
pl.save("figures/percent_growing.svg", width=5, height=5)

# Show some representative correlations
comp = both[(both.tradeoff == "0.5")]
w = within_samples[within_samples.tradeoff == "0.5"]
w.index = w.id
comp = comp[comp.id.isin(w.id[(w.n >= 10)])]
comp["rho"] = w.loc[comp.id, "rho"].round(2).values
comp.id = comp.id + " (r=" + comp.rho.astype(str) + ")"

pl = (
    ggplot(comp, aes(x="rate", y="growth_rate"))
    + geom_point()
    + facet_wrap("~ id", nrow=3)
    + scale_x_log10()
    + scale_y_log10()
    + labs(x="replication rate [a.u.]", y="predicted growth rate [1/h]")
)
pl.save("figures/corr_examples.png", width=12, height=6, dpi=300)
