import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from adjustText import adjust_text
from scipy.stats import ttest_ind

sample_keep = ["run_accession", "subset", "status", "type"]

SCFAs = {"butyrate": "EX_but_", "acetate": "EX_ac_",
         "propionate": "EX_ppa_"}


def test_vs_ref(x, y, ref="CTRL"):
    groups = x[x != ref].unique()
    tests = [pd.Series(ttest_ind(y[x == ref], y[x == g],
                                 equal_var=False)) for g in groups]
    tests = pd.DataFrame(tests)
    tests.columns = ["statistic", "pvalue"]
    tests["group"] = groups
    return tests


def box_jitter(x, y, **kwargs):
    sns.boxplot(x=x, y=y, color="white")
    sns.stripplot(x=x, y=y, color="black")
    print(test_vs_ref(x, y))


def export_rates_plot(fluxes, groups, samples):
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
    grid = sns.FacetGrid(fluxes, col="subset", row="metabolite",
                         sharey=False, sharex=False)
    g = grid.map(box_jitter, "name", "flux", color="white")
    return g


media = pd.read_csv("../results/minimal_media.csv", index_col=0).fillna(0.0)
media["sample"] = media.index
media = media.melt(id_vars="sample", var_name="reaction", value_name="flux")
metabolites = pd.read_csv("../results/metabolites.csv", index_col=0)
media["id"] = media.reaction.str.lstrip("EX_")
media = pd.merge(media, metabolites, on="id")
samples = pd.read_csv("../recent.csv")[sample_keep]
samples = samples.rename(columns={"run_accession": "sample"})
media = pd.merge(media, samples, on="sample")

mat = media.pivot("id", "sample", "flux")
mat = mat.apply(lambda x: x / x.abs().max(), axis=1)
g = sns.clustermap(mat, cmap="RdBu", figsize=(40, 42))
g.ax_heatmap.set_xlabel("")
g.ax_heatmap.set_ylabel("")
plt.savefig("media.png")
plt.close()

fluxes = pd.read_csv("../results/min_media_fluxes.csv.gz", compression="gzip")
fluxes = fluxes.melt(id_vars=["sample", "compartment"], var_name="reaction",
                     value_name="flux")
fluxes = fluxes[fluxes.reaction.str.startswith("EX_") &
                (fluxes.compartment != "medium")].dropna()
fluxes["taxa"] = fluxes.compartment + "_" + fluxes["sample"]

samples = pd.read_csv("../recent.csv")[
    ["run_accession", "status", "subset", "type"]]
samples = samples.rename(columns={"run_accession": "sample"})
samples.index = samples["sample"]
genera = pd.read_csv("../genera.csv")[["samples", "name", "reads"]]
totals = genera.groupby("samples").reads.sum()
genera["relative"] = genera.reads / totals[genera.samples].values
fluxes = pd.merge(fluxes, genera, left_on=["sample", "compartment"],
                  right_on=["samples", "name"])
fluxes["tot_flux"] = fluxes.flux * fluxes.relative
plt.tight_layout()
print("Production rates:")
export_rates_plot(fluxes[fluxes.tot_flux < 0], SCFAs, samples)
plt.savefig("scfas_prod.svg")
plt.close()

print("Consumption rates:")
export_rates_plot(fluxes[fluxes.tot_flux > 0], SCFAs, samples)
plt.savefig("scfas_consumption.svg")
plt.close()

print("Net rates:")
export_rates_plot(fluxes, SCFAs, samples)
plt.savefig("scfas_net.svg")
plt.close()

scfa = []
for name, filt in SCFAs.items():
    fl = fluxes[fluxes.reaction.str.contains(filt)].copy()
    fl = fl.groupby(["sample", "name"]).tot_flux.sum().reset_index()
    fl["metabolite"] = name
    scfa.append(fl)
    mat = fl.pivot("sample", "name", "tot_flux")
    mat = mat.loc[:, mat.abs().mean() > 0.5].fillna(0)
    sns.clustermap(mat, cmap="seismic", center=0, figsize=(mat.shape[1]/2, 10),
                   yticklabels=False)
    plt.savefig(name + ".svg")
    plt.close()
scfa = pd.concat(scfa)
ord = scfa.groupby(["name"]).tot_flux.apply(lambda x: x.abs().mean())
ord = ord.sort_values(ascending=False)
plt.figure(figsize=(8, 7))
plt.axhline(0, c="black")
ax = sns.pointplot(x="name", y="tot_flux", hue="metabolite",
                   data=scfa[scfa.name.isin(ord.index[ord > 1])], ci="sd",
                   order=ord.index[ord > 1], join=False, dodge=True)
ax.grid(axis="x", color="gainsboro")
ax.xaxis.set_tick_params(rotation=90)
sns.despine(left=True, bottom=True)
plt.xlabel("")
plt.tight_layout()
plt.savefig("scfas.svg")
plt.close()


mat = fluxes.pivot("taxa", "reaction", "flux").fillna(0.0)
taxa = mat.index.str.split("_").str[0]
tsne = TSNE(n_components=2).fit_transform(mat)
tsne = pd.DataFrame(tsne, columns=["x", "y"], index=mat.index)
tsne["taxa"] = taxa
sns.set(font_scale=1.5, style="ticks")
g = sns.FacetGrid(tsne, hue="taxa", size=10, aspect=16/10)
gm = g.map(plt.scatter, "x", "y", alpha=0.25)
means = tsne.groupby(taxa).agg("median").reset_index()
texts = means.apply(lambda df: plt.text(df.x, df.y, df.taxa, alpha=0.65),
                    axis=1)
texts = adjust_text(texts, force_text=(0.02, 0.5),
                    arrowprops=dict(arrowstyle='-|>', alpha=0.5, color="k"))
plt.savefig("individual_media.png", dpi=200)
plt.close()
