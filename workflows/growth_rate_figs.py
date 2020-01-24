import pandas as pd
import numpy as np
from plotnine import *
from scipy.stats import pearsonr

theme_set(theme_minimal())

rates = (
    pd.read_csv("data/tradeoff.csv")
    .dropna(subset=["abundance"])
    .rename(columns={"sample": "id", "compartments": "genus"})
)

rates["log_rates"] = np.log10(rates.growth_rate)
rates.loc[rates.growth_rate <= 0, "log_rates"] = -16

pos = rates.query("growth_rate > 1e-6 and tradeoff == 0.5")
o = pos.groupby("genus").log_rates.mean().sort_values().index
pos.genus = pd.Categorical(pos.genus, categories=o, ordered=True)

pl = (
    ggplot(pos, aes(x="genus", y="growth_rate", color="genus"))
    + geom_jitter(width=0.2, height=0, alpha=0.5, stroke=0, size=2)
    + stat_summary(fun_y=np.mean, geom="point", size=2, fill="white", stroke=1)
    + scale_y_log10()
    + guides(color=False)
    + labs(x="", y="growth rate [1/h]")
    + theme(axis_text_x=element_text(rotation=45, hjust=1, vjust=1))
)
pl.save("figures/gcs.svg", width=12, height=4)

community_growth = (
    rates[rates.tradeoff == 0.5]
    .groupby("id")
    .apply(lambda df: sum(df.abundance * df.growth_rate))
)
print(community_growth.describe())


rs = rates[rates.tradeoff == 0.5]
scale = rs.groupby("id").abundance.apply(
    lambda a: community_growth[a.name] / a.pow(2).sum()
)
x = np.linspace(0, 1, 100)
free = pd.DataFrame(
    {
        "abundance": x,
        "max_rate": scale.max() * x,
        "min_rate": scale.min() * x,
        "mean_rate": scale.mean() * x,
    }
)
pl = (
    ggplot(rs, aes(x="abundance", y="growth_rate"))
    + geom_point(stroke=0, alpha=0.2)
    + labs(x="", y="")
    + theme_classic()
)
pl.save("figures/rate_vs_abundance.png", width=3, height=3, dpi=300)

subset = rs.id.unique()[0:10]
free_subset = (
    rs[rs.id.isin(subset)]
    .groupby("id")
    .apply(
        lambda df: pd.DataFrame(
            {"abundance": x, "growth_rate": scale[df.name] * x}
        )
    )
    .reset_index()
)
pl = (
    ggplot(
        rs[rs.id.isin(subset)], aes(x="abundance", y="growth_rate", color="id")
    )
    + geom_point(size=2)
    + geom_line(data=free_subset, linetype="dashed")
    + labs(x="abundance [%]", y="growth rate [1/h]")
    + guides(color=False)
    + theme_classic()
)
pl.save("figures/rate_vs_abundance_individual.svg", width=6, height=5)

