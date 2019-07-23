import pandas as pd
import micom

keep = [
    "samples",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
    "oxygenstat",
    "metabolism",
    "gram",
    "mtype",
    "genes",
    "file",
    "reads",
    "relative",
]


def reduce_group(df):
    new = df.iloc[0, :]
    new["file"] = "|".join(df.id.apply(lambda id: f"{id}.xml"))
    return new


agora = micom.data.agora
agora = agora.rename(columns={"mclass": "class", "superkingdom": "kingdom"})
agora = agora.groupby("species").apply(reduce_group).reset_index(drop=True)


species = pd.read_csv("data/abundances.csv")
species.species = species.species.str.split(" ").str[1]

species_models = pd.merge(
    species,
    agora,
    on=["kingdom", "phylum", "class", "order", "family", "genus", "species"],
)
species_models = species_models.rename(columns={"id_x": "samples"})[keep]
species_models.to_csv("data/species.csv", index=False)
