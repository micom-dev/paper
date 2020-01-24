"""Calculate stats for taxa assignments and availability in AGORA."""

import pandas as pd
import micom


tax = pd.read_csv("data/abundances.csv").query("kingdom == 'Bacteria'")
tax.relative = tax.groupby("id").relative.apply(lambda a: a / a.sum())
tax.taxa_id = tax.taxa_id.str.replace("*", "").astype("int")
agora = micom.data.agora
agora.species = agora.genus + " " + agora.species
taxa = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]


def taxa_stats(taxonomy, rank, agora):
    res = pd.Series(
        index=[
            "n_unique",
            "mean_percent_assigned",
            "sd_percent_assigned",
            "n_model",
            "mean_percent_model",
            "sd_percent_model",
        ]
    )
    assigned = taxonomy.dropna(subset=[rank]).groupby("id").relative.sum()
    if rank == "superkingdom":
        arank = "kingdom"
    elif rank == "class":
        arank = "mclass"
    else:
        arank = rank
    good = taxonomy[rank].isin(agora[arank].dropna())
    has_model = taxonomy[good].groupby("id").relative.sum()
    res.iloc[:] = [
        taxonomy[rank].nunique(),
        assigned.mean(),
        assigned.std(),
        taxonomy[good][rank].nunique(),
        has_model.mean(),
        has_model.std(),
    ]
    res.name = rank
    return res


stats = pd.concat([taxa_stats(tax, ta, agora) for ta in taxa], axis=1)
print(stats)
