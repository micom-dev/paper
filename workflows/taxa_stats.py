"""Calculate stats for taxa assignments and availability in AGORA."""

import pandas as pd


tax = pd.read_csv("results/taxonomy.csv")
agora = pd.read_csv("agora.csv")
taxa = tax["rank"].unique()

def taxa_stats(taxonomy):
    res = pd.Series(index=["n_unique",
                           "mean_percent_assigned",
                           "sd_percent_assigned",
                           "mean_percent_model",
                           "sd_percent_model"])
    assigned = taxonomy.dropna().groupby("id").relative.sum()
    rank = taxonomy["rank"].iloc[0]
    if rank in ["genus", "species"]:
        good = taxonomy.taxid.isin(agora["ncbi_" + rank + "_id"])
    elif rank == "superkingdom":
        good = taxonomy.name.isin(agora["kingdom"])
    else:
        good = taxonomy.name.isin(agora[rank])
    has_model = taxonomy[good].dropna().groupby("id").relative.sum()
    res.iloc[:] = [taxonomy.name.nunique(),
                  assigned.mean(), assigned.std(),
                  has_model.mean(), has_model.std()]
    res.name = rank
    return res


stats = pd.concat([taxa_stats(tax[tax["rank"] == ta]) for ta in taxa],
                  axis=1)

