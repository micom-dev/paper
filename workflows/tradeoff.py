"""Get growth rates over a variety of tradeoffs."""

import micom
from micom import load_pickle
from micom.workflows import workflow
import numpy as np
import pandas as pd


tradeoffs = np.arange(0.05, 1.01, 0.05)
micom.logger.file_logger("micom.log")
logger = micom.logger.logger
try:
    max_procs = snakemake.threads
except NameError:
    max_procs = 20
logger.info("Using %d threads..." % max_procs)


def growth_rates(sam):
    com = load_pickle("data/models/" + sam + ".pickle")
    sol = com.optimize()
    rates = sol.members
    rates["tradeoff"] = np.nan
    rates["sample"] = sam
    df = [rates]

    # Get growth rates
    try:
        sol = com.cooperative_tradeoff(fraction=tradeoffs)
    except Exception as e:
        logger.warning("Sample %s could not be optimized\n %s" % (sam, str(e)))
        return pd.DataFrame(
            {"tradeoff": tradeoffs, "growth_rate": np.nan, "sample": sam}
        )
    for i, s in enumerate(sol.solution):
        rates = s.members
        rates["tradeoff"] = tradeoffs[i]
        rates["sample"] = sam
        df.append(rates)
    df = pd.concat(df)
    return df


samples = pd.read_csv("data/recent.csv")
rates = workflow(growth_rates, samples.run_accession, max_procs)
rates = pd.concat(rates)
rates.to_csv("data/tradeoff.csv")
