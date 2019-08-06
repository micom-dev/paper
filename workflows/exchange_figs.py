import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from adjustText import adjust_text
from scipy.stats import ttest_ind
import numpy as np

sample_keep = ["run_accession", "subset", "status", "type"]

SCFAs = {"butyrate": "EX_but_", "acetate": "EX_ac_", "propionate": "EX_ppa_"}


def test_vs_ref(x, y, ref="CTRL"):
    groups = x[x != ref].unique()
    tests = [
        pd.Series(ttest_ind(y[x == ref], y[x == g], equal_var=False))
        for g in groups
    ]
    tests = pd.DataFrame(tests)
    tests.columns = ["statistic", "pvalue"]
    tests["group"] = groups
    return tests


def box_jitter(x, y, **kwargs):
    sns.boxplot(x=x, y=y, color="white")
    sns.stripplot(x=x, y=y, color="black")
    print(test_vs_ref(y, x))


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
    grid = sns.FacetGrid(
        fluxes,
        col="subset",
        row="metabolite",
        sharey=False,
        sharex=True,
        aspect=1.6,
    )
    g = grid.map(box_jitter, "flux", "name", color="white")
    return g


media = pd.read_csv("data/minimal_imports.csv", index_col=0).fillna(0.0)
media["sample"] = media.index
media = media.melt(id_vars="sample", var_name="reaction", value_name="flux")
metabolites = pd.read_csv("data/metabolites.csv")
metabolites["id"] = metabolites.abbreviation + "_m"
media["id"] = media.reaction.str.lstrip("EX_")
media = pd.merge(media, metabolites, on="id")
samples = pd.read_csv("data/recent.csv")[sample_keep]
samples = samples.rename(columns={"run_accession": "sample"})
media = pd.merge(media, samples, on="sample")

mat = media.pivot("id", "sample", "flux")
mat = mat.apply(lambda x: x / x.abs().max(), axis=1)
g = sns.clustermap(mat, cmap="RdBu", figsize=(60, 70), center=0)
g.ax_heatmap.set_xlabel("")
g.ax_heatmap.set_ylabel("")
plt.savefig("figures/media.png")
plt.close()

fluxes = pd.read_csv("data/minimal_fluxes.csv.gz", compression="gzip")
fluxes = fluxes.melt(
    id_vars=["sample", "compartment"], var_name="reaction", value_name="flux"
)
fluxes = fluxes[
    fluxes.reaction.str.startswith("EX_") & (fluxes.compartment != "medium")
].dropna()
fluxes["taxa"] = fluxes.compartment + "_" + fluxes["sample"]
fluxes["name"] = fluxes.compartment.str.replace("_", " ")

samples = pd.read_csv("data/recent.csv")[
    ["run_accession", "status", "subset", "type"]
]
samples = samples.rename(columns={"run_accession": "sample"})
samples.index = samples["sample"]
species = pd.read_csv("data/genera.csv")[["samples", "genus", "reads"]]
species["name"] = species.genus
totals = species.groupby("samples").reads.sum()
species["relative"] = species.reads / totals[species.samples].values
fluxes = pd.merge(
    fluxes, species, left_on=["sample", "name"], right_on=["samples", "name"]
)
fluxes["tot_flux"] = fluxes.flux * fluxes.relative
plt.tight_layout()
print("Production rates:")
export_rates_plot(fluxes[fluxes.tot_flux < 0], SCFAs, samples)
plt.savefig("figures/scfas_prod.svg")
plt.close()

print("Consumption rates:")
export_rates_plot(fluxes[fluxes.tot_flux > 0], SCFAs, samples)
plt.savefig("figures/scfas_consumption.svg")
plt.close()

print("Net rates:")
export_rates_plot(fluxes, SCFAs, samples, log=True)
plt.savefig("figures/scfas_net.svg")
plt.close()

scfa = []
for name, filt in SCFAs.items():
    fl = fluxes[fluxes.reaction.str.contains(filt)].copy()
    fl = fl.groupby(["sample", "name"]).tot_flux.sum().reset_index()
    fl["metabolite"] = name
    scfa.append(fl)
    mat = fl.pivot("sample", "name", "tot_flux")
    mat = mat.loc[
        :, mat.abs().mean().sort_values(ascending=False).index[0:5]
    ].fillna(0)
    sns.clustermap(
        mat, cmap="seismic", center=0, figsize=(6, 12), yticklabels=False
    )
    plt.savefig("figures/" + name + ".svg")
    plt.close()
scfa = pd.concat(scfa)
ord = scfa.groupby(["name"]).tot_flux.apply(lambda x: x.abs().mean())
ord = ord.sort_values(ascending=False)
plt.figure(figsize=(8, 7))
plt.axhline(0, c="black")
ax = sns.pointplot(
    x="name",
    y="tot_flux",
    hue="metabolite",
    data=scfa[scfa.name.isin(ord.index[ord > 1])],
    ci="sd",
    order=ord.index[ord > 1],
    join=False,
    dodge=True,
)
ax.grid(axis="x", color="gainsboro")
ax.xaxis.set_tick_params(rotation=90)
sns.despine(left=True, bottom=True)
plt.xlabel("")
plt.tight_layout()
plt.savefig("figures/scfas.svg")
plt.close()


mat = fluxes.pivot("taxa", "reaction", "flux").fillna(0.0)
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
