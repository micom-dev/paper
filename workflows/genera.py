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
    "oxygen_status",
    "metabolism",
    "gram",
    "type",
    "pubseed_id",
    "genome_size",
    "genes",
    "file",
    "name",
    "reads",
    "relative",
]


def reduce_group(df):
    new = df.iloc[0, :]
    new["file"] = "|".join(df.id.apply(lambda id: f"{id}.xml"))
    return new


agora = micom.data.agora
agora_genus = agora.groupby("genus").apply(reduce_group)

genera = pd.read_csv("data/abundances.csv")
genera = genera[genera["rank"] == "genus"]

genus_models = pd.merge(genera, agora_genus, on="genus")
genus_models = genus_models.rename(columns={"id_x": "samples"})[keep]
genus_models.to_csv("data/genera.csv", index=False)
