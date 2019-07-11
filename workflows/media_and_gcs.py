"""Extract minimal growth media and growth rates."""

import pandas as pd
from micom import load_pickle
from micom.media import minimal_medium
from micom.workflows import workflow


max_procs = 20


def media_and_gcs(sam):
    com = load_pickle("models/" + sam + ".pickle")

    # Get growth rates
    sol = com.cooperative_tradeoff(fraction=0.9)
    rates = sol.members["growth_rate"].copy()
    rates["community"] = sol.growth_rate
    rates.name = s

    # Get the minimal medium
    med = minimal_medium(com, 0.95 * sol.growth_rate)
    med.name = s
    return {"medium": med, "gcs": rates}


samples = pd.read_csv("recent.csv")
gcs = pd.DataFrame()
media = pd.DataFrame()

results = workflow(media_and_gcs, samples.run_accession, max_procs)

for s in results:
    gcs = gcs.append(results["gcs"])
    media = media.append(results["media"])

gcs.to_csv("data/growth_rates.csv")
media.to_csv("data/minimal_media.csv")
