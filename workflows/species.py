import pandas as pd
import micom

keep = [
    "samples",
    "kingdom_y",
    "phylum_y",
    "class_y",
    "order_y",
    "family_y",
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
    "taxa_id_y",
]


def reduce_group(df):
    new = df.iloc[0, :]
    new["file"] = "|".join(df.id.apply(lambda id: f"{id}.xml"))
    return new


agora = micom.data.agora
agora = agora.rename(columns={"mclass": "class", "superkingdom": "kingdom"})
agora = (
    agora.groupby(["genus", "species"])
    .apply(reduce_group)
    .reset_index(drop=True)
)


species = pd.read_csv("data/abundances.csv")
species.species = species.species.str.split(" ").str[1]
species.taxa_id = species.taxa_id.str.replace("*", "").astype("int")

species_models = pd.merge(species, agora, on=["genus", "species"])
species_models = species_models.rename(columns={"id_x": "samples"})[keep]
species_models.columns = species_models.columns.str.replace("_y", "")
species_models.to_csv("data/species.csv", index=False)
