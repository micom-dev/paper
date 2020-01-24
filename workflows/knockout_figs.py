"""Analyze interactions between taxa."""

import pandas as pd
import networkx as nx
import nxviz
import matplotlib.pyplot as plt
from plotnine import *

theme_set(theme_minimal())

ko = pd.read_csv("data/knockouts.csv", index_col=0)
ko["knocked"] = ko.index
ko = ko.melt(
    id_vars=["knocked", "sample"], var_name="genus", value_name="change"
).dropna()

pos = ko.groupby(["knocked", "genus"]).change.mean().reset_index()
pos = pos[pos.knocked != pos.genus].sort_values(by="knocked")
counts = (
    ko[ko.change.abs() > 1e-2]
    .groupby("knocked")
    .apply(
        lambda df: pd.DataFrame(
            {
                "counts": [
                    df.change[df.change > 0].count(),
                    df.change[df.change < 0].count(),
                ],
                "type": ["competitive", "cooperative"],
            }
        )
    )
).reset_index()
counts = counts.sort_values(by="counts", ascending=False)
counts.knocked = pd.Categorical(
    counts.knocked, counts.knocked.unique()[::-1], ordered=True
)
print(counts.groupby("type").sum())

pl = (
    ggplot(counts, aes(x="knocked", y="counts", fill="type"))
    + geom_col()
    + labs(x="", y="strong interactions")
    + coord_flip()
)
pl.save("figures/knockout_counts.svg", width=4, height=12)

graph = nx.from_pandas_edgelist(
    pos[pos.change.abs() > 1e-2], "knocked", "genus", "change"
)

for idx, _ in graph.nodes(data=True):
    graph.node[idx]["degree"] = float(graph.degree(idx))
circos = nxviz.CircosPlot(
    graph,
    node_labels=True,
    rotate_labels=True,
    edge_cmap="bwr",
    edge_color="change",
    edge_limits=(-pos.change.abs().max(), pos.change.abs().max()),
    node_color="degree",
    figsize=(9, 7.5),
)
# circos.figure.set_dpi(200)
circos.draw()
plt.tight_layout(rect=[0.05, 0.15, 0.7, 0.8])
plt.savefig("figures/circos.svg")
plt.close()

ns = (
    ko[(ko.change.abs() > 1e-2) & (ko.knocked != ko.genus)]
    .groupby(["genus", "knocked"])
    .sample.count()
    .reset_index()
)

pl = (
    ggplot(ns, aes(x="sample"))
    + geom_histogram(bins=32, fill="lightgray", color="black")
    + labs(x="observed samples")
)
pl.save("figures/interaction_prevalence.svg", width=5, height=5)
