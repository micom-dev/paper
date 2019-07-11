"""Extract minimal growth media and growth rates."""

import pandas as pd
import micom
from micom import load_pickle
from micom.media import minimal_medium
from micom.workflows import workflow


micom.logger.file_logger("micom.log")
logger = micom.logger.logger
try:
    max_procs = snakemake.threads
except NameError:
    max_procs = 20


def media_and_gcs(sam):
    com = load_pickle("models/" + sam + ".pickle")

    # Get growth rates
    sol = com.cooperative_tradeoff(fraction=0.5)
    rates = sol.members["growth_rate"].copy()
    rates["community"] = sol.growth_rate
    rates.name = sam

    # Get the minimal medium
    med = minimal_medium(com, 0.95 * sol.growth_rate)
    med.name = sam

    # Apply medium and reoptimize
    com.medium = med
    sol = com.cooperative_tradeoff(fraction=0.5, fluxes=True, pfba=False)
    fluxes = sol.fluxes
    fluxes["sample"] = sam
    return {"medium": med, "gcs": rates, "fluxes": fluxes}


samples = pd.read_csv("recent.csv")
gcs = pd.DataFrame()
media = pd.DataFrame()
fluxes = pd.DataFrame()

results = workflow(media_and_gcs, samples.run_accession, max_procs)

for r in results:
    gcs = gcs.append(r["gcs"])
    media = media.append(r["media"])
    fluxes = fluxes.append(r["fluxes"])

gcs.to_csv("data/growth_rates.csv")
media.to_csv("data/minimal_imports.csv")
fluxes.to_csv("data/minimal_fluxes.csv.gz", compression="gzip")
