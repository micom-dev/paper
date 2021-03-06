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
agora_genus = agora.groupby("genus").apply(reduce_group).reset_index(drop=True)

genera = pd.read_csv("data/abundances.csv")
genera = genera.groupby(["id", "genus"]).sum().reset_index()

genus_models = pd.merge(genera, agora_genus, on="genus", suffixes=["_x", ""])
genus_models = genus_models.rename(
    columns={"id_x": "samples", "mclass": "class"}
)[keep]
genus_models.to_csv("data/genera.csv", index=False)
