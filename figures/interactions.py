"""Analyze interactions between taxa."""

import pandas as pd
import networkx as nx
import nxviz
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

ko = pd.read_csv("../knockouts.csv", index_col=0)
ko["knocked"] = ko.index
ko = ko.melt(id_vars=["knocked", "sample"], var_name="genus",
             value_name="change").dropna()
maxs = ko.groupby("genus").change.apply(lambda x: x.abs().max())
normalized = ko.change.values / maxs[ko.genus]
normalized.index = ko.index
ko["norm_change"] = normalized


pos = ko[ko.norm_change.abs() > 0.5]
pos = pos[pos.knocked != pos.genus]
graph = nx.from_pandas_dataframe(pos,
                                 "knocked", "genus", "norm_change")

for idx, _ in graph.nodes(data=True):
    graph.node[idx]['degree'] = float(graph.degree(idx))
circos = nxviz.CircosPlot(graph, node_labels=True, rotate_labels=True,
                          edge_cmap="bwr",
                          edge_color="norm_change", node_color="degree",
                          figsize=(9, 7.5))
#circos.figure.set_dpi(200)
circos.draw()
plt.tight_layout(rect=[0.05, 0.15, 0.7, 0.8])
plt.savefig("circos.svg")
