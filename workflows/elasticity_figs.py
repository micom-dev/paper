"""Analyze elasticities for reference samples."""

import pandas as pd
import networkx as nx
import nxviz
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

samples = {
    "ERR260275": "healthy",
    "ERR260214": "T2D metformin+",
    "ERR260174": "T2D metformin-",
}
metabolites = pd.read_csv("data/metabolites.csv")[
    ["abbreviation", "fullName", "keggId"]
]
SCFAs = {"butyrate": "EX_but(e)", "acetate": "EX_ac(e)", "propionate": "EX_ppa(e)"}


def direction(els, rid):
    return els[els.reaction == rid].direction.unique()[0]


elast = []
for sa in samples:
    e = pd.read_csv("data/elasticities_" + sa + ".csv")
    e["id"] = sa
    elast.append(e)
elast = pd.concat(elast)
elast = elast[elast.direction == "forward"]
elast["scfa"] = float("nan")
for name, pattern in SCFAs.items():
    elast.loc[elast.reaction.str.startswith(pattern), "scfa"] = name

production = (
    elast.groupby(["id", "effector", "scfa"]).elasticity.sum().reset_index()
)
production["type"] = production.effector.apply(
    lambda e: "diet" if e.startswith("EX_") else "abundance"
)
production["scfa"] = production.scfa + " " + production.id.map(samples)
production["abbreviation"] = production.effector.str.replace("(EX_)|(_m)", "")
production = pd.merge(production, metabolites, on="abbreviation", how="left")
production.loc[production.fullName.isna(), "fullName"] = production.loc[
    production.fullName.isna(), "abbreviation"
]

cmap = sns.color_palette()[0:2]
cmap = dict(zip(production.type.unique(), cmap))
production.fullName = production.fullName.apply(
    lambda s: s if len(s) < 24 else s[:24] + "..."
)
ty = production[["fullName", "type"]].drop_duplicates()
ty = pd.Series(ty.type.values, index=ty.fullName)
typecols = ty.map(cmap).rename("type")
production = production.sort_values(by="id")
mat = production.pivot_table(
    index="scfa", columns="fullName", values="elasticity", fill_value=0
)
g = sns.clustermap(
    mat,
    cmap="seismic",
    vmin=-production.elasticity.abs().max(),
    vmax=production.elasticity.abs().max(),
    figsize=(32, 4),
    col_colors=typecols,
    yticklabels=True,
    xticklabels=True,
    row_cluster=False,
    cbar_kws={"fraction": 1.0},
)
for label in cmap:
    g.ax_col_dendrogram.bar(0, 0, color=cmap[label], label=label, linewidth=0)
g.ax_col_dendrogram.legend(loc="best", ncol=3, frameon=False)
g.ax_heatmap.set_xlabel("")
g.ax_heatmap.set_ylabel("")
cold = g.ax_col_dendrogram.get_position()
g.ax_col_dendrogram.set_position(
    [cold.x0, cold.y0, cold.width, cold.height * 5]
)
g.ax_row_dendrogram.set_visible(False)
g.cax.set_position([0.25, 0.5, 0.01, 0.5])
g.savefig("figures/elasticities.png", dpi=300)

but = production[production.scfa.str.startswith("butyrate")].sort_values(
    by="elasticity", ascending=False
)
print(but.head(20))
