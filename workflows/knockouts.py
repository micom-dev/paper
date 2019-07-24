"""Perform taxa knockouts."""

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


def knockout(sam):
    com = load_pickle("data/models/" + sam + ".pickle")
    ko = com.knockout_species(fraction=0.5)
    ko["sample"] = sam
    return ko


samples = pd.read_csv("data/recent.csv")
kos = workflow(knockout, samples.run_accession, max_procs)
pd.concat(kos).to_csv("data/knockouts.csv")
