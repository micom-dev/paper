import pandas as pd
from plotnine import *
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from adjustText import adjust_text
from scipy.stats import ttest_ind
import numpy as np
from itertools import combinations
from scipy.spatial.distance import pdist
from skbio.stats.distance import DistanceMatrix, permanova

sample_keep = ["run_accession", "subset", "status", "type"]
SCFAs = {
    "butyrate": "EX_but\\(e\\)",
    "acetate": "EX_ac\\(e\\)",
    "propionate": "EX_ppa\\(e\\)",
}
theme_set(theme_minimal())


def pairwise_tests(x, y, name="tests", ref="CTRL"):
    groups = x.unique()
    combs = combinations(groups, 2)
    tests = [
        pd.Series(
            ttest_ind(
                y[x == g[0]].dropna(), y[x == g[1]].dropna(), equal_var=False
            )
            + (len(y), g[0], g[1]),
            index=["t", "pval", "n", "x", "y"],
        )
        for g in combs
    ]
    tests = pd.DataFrame(tests)
    tests["name"] = name
    return tests


def export_rates_plot(fluxes, groups, samples, log=False):
    dfs = []
    for name, filt in groups.items():
        df = fluxes[fluxes.reaction.str.contains(filt)].copy()
        res = samples.copy()
        df = df.groupby(["sample", "compartment"]).tot_flux.sum().reset_index()
        res["flux"] = df.groupby("sample").tot_flux.sum().abs()
        res["metabolite"] = name
        dfs.append(res)
    fluxes = pd.concat(dfs)
    fluxes.loc[fluxes.status == "ND", "status"] = ""
    fluxes["name"] = fluxes.status + " " + fluxes.type.fillna("")
    fluxes.name = fluxes.name.str.strip()
    fluxes = fluxes.sort_values("name")
    if log:
        fluxes.flux = np.log10(fluxes.flux)
    fluxes.groupby(["metabolite", "subset"]).apply(
        lambda df: print(
            pairwise_tests(df["name"], df.flux, df.name[0] + " " + df.name[1])
        )
    )
    pl = (
        ggplot(fluxes, aes(x="name", y="flux", color="metabolite"))
        + geom_boxplot(outlier_color="none")
        + geom_jitter(width=0.2, height=0)
        + facet_grid("metabolite ~ subset", scales="free")
        + guides(color=False)
        + labs(x="", y="flux [mmol/h]")
        + theme(axis_text_x=element_text(rotation=45, vjust=1, hjust=1))
    )
    return pl


samples = pd.read_csv("data/recent.csv")[
    ["run_accession", "status", "subset", "type"]
]
samples = samples.rename(columns={"run_accession": "sample"})
samples.index = samples["sample"]

media = pd.read_csv("data/minimal_imports.csv", index_col=0).fillna(0.0)
media["sample"] = media.index
media = media.melt(id_vars="sample", var_name="reaction", value_name="flux")
media = media[media.flux > 0]
metabolites = pd.read_csv("data/metabolites.csv")
metabolites["id"] = metabolites.abbreviation + "_m"
media["id"] = media.reaction.str.lstrip("EX_")
media = pd.merge(media, metabolites, on="id")

mat = media.dropna(subset=["flux"]).pivot("id", "sample", "flux").fillna(0)
stype = (
    samples.loc[mat.columns].status.fillna("")
    + " "
    + samples.loc[mat.columns].type.fillna("")
).replace({"ND CTRL": "CTRL"})
m = dict(zip(stype.unique(), sns.color_palette()[0 : (stype.nunique() + 1)]))
stype_cols = stype.map(m)
stype_cols.name = "group"
mat = mat.apply(lambda x: x / x.abs().max(), axis=1)
g = sns.clustermap(
    mat,
    cmap="BuPu",
    figsize=(18, 18),
    col_colors=stype_cols,
    yticklabels=True,
    xticklabels=True,
)
for label in m:
    g.ax_col_dendrogram.bar(0, 0, color=m[label], label=label, linewidth=0)
g.ax_col_dendrogram.legend(loc="best", ncol=2, frameon=False)
g.ax_heatmap.set_xlabel("")
g.ax_heatmap.set_ylabel("")
col = g.ax_col_colors.get_position()
g.ax_col_colors.set_position([col.x0, col.y0, col.width, col.height * 0.3])
cold = g.ax_col_dendrogram.get_position()
g.ax_col_dendrogram.set_position(
    [cold.x0, cold.y0 - col.height * 0.7, cold.width, cold.height]
)
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize=6)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize=6)
g.ax_col_colors.set_yticklabels(
    g.ax_col_colors.get_ymajorticklabels(), fontsize=12
)
g.savefig("figures/media.png", dpi=300)
plt.close()

fluxes = pd.read_csv("data/minimal_fluxes.csv.gz", compression="gzip")
fluxes = fluxes.melt(
    id_vars=["sample", "compartment"], var_name="reaction", value_name="flux"
)
fluxes = fluxes[
    fluxes.reaction.str.startswith("EX_")
    & (fluxes.compartment != "medium")
    & fluxes.reaction.str.endswith("(e)")
].fillna(0)
fluxes["taxa"] = fluxes.compartment + "_" + fluxes["sample"]
fluxes["name"] = fluxes.compartment.str.replace("_", " ")

species = pd.read_csv("data/genera.csv")[["samples", "genus", "reads"]]
species["name"] = species.genus
totals = species.groupby("samples").reads.sum()
species["relative"] = species.reads / totals[species.samples].values
fluxes = pd.merge(
    fluxes, species, left_on=["sample", "name"], right_on=["samples", "name"]
)
fluxes["tot_flux"] = fluxes.flux * fluxes.relative
print("Production rates:")
pl = export_rates_plot(fluxes[fluxes.tot_flux > 0], SCFAs, samples)
pl.save("figures/scfas_prod.svg", width=4, height=6)

print("Consumption rates:")
pl = export_rates_plot(fluxes[fluxes.tot_flux < 0], SCFAs, samples)
pl.save("figures/scfas_consumption.svg", width=4, height=6)

print("Net rates:")
pl = export_rates_plot(fluxes, SCFAs, samples)
pl.save("figures/scfas_net.svg", width=4, height=6)

scfa = []
for name, filt in SCFAs.items():
    fl = fluxes[fluxes.reaction.str.contains(filt)].copy()
    fl = fl.groupby(["sample", "name", "relative"]).flux.sum().reset_index()
    fl["metabolite"] = name
    scfa.append(fl)
    mat = fl.pivot("sample", "name", "flux")
    mat = mat.loc[
        :, mat.abs().mean().sort_values(ascending=False).index[0:5]
    ].fillna(0)
    sns.clustermap(
        mat, cmap="seismic", center=0, figsize=(6, 12), yticklabels=False
    )
    plt.savefig("figures/" + name + ".svg")
    plt.close()
scfa = pd.concat(scfa)
ord = scfa.groupby(["name"]).relative.mean().sort_values(ascending=False)
scfa.name = pd.Categorical(
    scfa.name.astype("str"), categories=ord.index[::-1], ordered=True
)
pl = (
    ggplot(
        scfa[scfa.name.isin(ord[ord > 0.01].index)],
        aes(x="name", y="flux", fill="metabolite", group="name"),
    )
    + geom_hline(yintercept=0, linetype="dashed")
    + geom_boxplot(outlier_color="none")
    + facet_wrap("~ metabolite", scales="free", nrow=3)
    + guides(fill=False)
    + labs(x="", y="flux [mmol/(h·gDW)]")
    + coord_flip()
)
pl.save("figures/scfas.svg", width=2.5, height=15)

pl = (
    ggplot(
        scfa[scfa.name.isin(ord[ord > 0.01].index)],
        aes(x="name", y="relative"),
    )
    + geom_boxplot(outlier_color="none")
    + labs(x="", y="relative\nabundance")
    + coord_flip()
)
pl.save("figures/relative_abundance.svg", width=2, height=5)


mat = fluxes[fluxes.flux < 0].pivot("taxa", "reaction", "flux").fillna(0.0)
taxa = mat.index.str.split("_ERR").str[0]
tsne = TSNE(n_components=2).fit_transform(mat)
tsne = pd.DataFrame(tsne, columns=["x", "y"], index=mat.index)
tsne["taxa"] = taxa
sns.set(font_scale=1.5, style="ticks")
g = sns.FacetGrid(tsne, hue="taxa", height=10, aspect=16 / 10)
gm = g.map(plt.scatter, "x", "y", alpha=0.25)
means = tsne.groupby(taxa).agg("median").reset_index()
texts = means.apply(
    lambda df: plt.text(df.x, df.y, df.taxa, alpha=0.65), axis=1
)
texts = adjust_text(
    texts,
    force_text=(0.02, 0.5),
    arrowprops=dict(arrowstyle="-|>", alpha=0.5, color="k"),
)
plt.savefig("figures/individual_media.png", dpi=200)
plt.close()

# Some statistics about metabolite usage
# indicator matrix 0 = metabolite not consumed, 1 = metabolite consumed
binary = mat.where(mat < -1e-6, 0).where(mat > -1e-6, 1)

# Jaccard distances = 1 - percent overlap
J = pdist(binary, "jaccard")
print("Jaccard distances:", pd.Series(J).describe(), sep="\n")

# euclidean distances
E = pdist(mat, "euclidean")

# Test whether genus explains a good amount of that variation
p = permanova(DistanceMatrix(E), taxa)
r2 = 1 - 1 / (1 + p[4] * p[3] / (p[2] - p[3] - 1))
p["R2"] = r2
print("PERMANOVA on euclidean distances:", p, sep="\n")
