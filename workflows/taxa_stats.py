"""Calculate stats for taxa assignments and availability in AGORA."""

import pandas as pd
import micom


tax = pd.read_csv("data/abundances.csv")
tax.species = tax.species.str.split(" ").str[1]
agora = micom.data.agora
taxa = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]


def taxa_stats(taxonomy, rank):
    res = pd.Series(
        index=[
            "n_unique",
            "mean_percent_assigned",
            "sd_percent_assigned",
            "mean_percent_model",
            "sd_percent_model",
        ]
    )
    assigned = taxonomy.dropna().groupby("id").relative.sum()
    if rank == "superkingdom":
        arank = "kingdom"
    if rank == "class":
        arank = "mclass"
    else:
        arank = rank
    good = taxonomy[rank].isin(agora[arank])
    has_model = taxonomy[good].dropna().groupby("id").relative.sum()
    res.iloc[:] = [
        taxonomy[rank].nunique(),
        assigned.mean(),
        assigned.std(),
        has_model.mean(),
        has_model.std(),
    ]
    res.name = rank
    return res


stats = pd.concat([taxa_stats(tax, ta) for ta in taxa], axis=1)

