"""Calculate the elasticities for a set of built models."""

from os.path import isfile
import micom
from micom import load_pickle
from micom.elasticity import elasticities
from micom.workflows import workflow


logger = micom.logger.logger
logger.add("micom.log")
try:
    max_procs = snakemake.threads
except NameError:
    max_procs = 20


def elasticity_worker(sam):
    """Get the exchange elasticities for a sample."""
    model_file = "data/models/" + sam + ".pickle"
    out_file = "data/elasticities_" + sam + ".csv"
    if isfile(out_file):
        return
    com = load_pickle(model_file)
    elast = elasticities(com, fraction=0.5)
    elast.to_csv(out_file, index=False)


samples = ["ERR260275", "ERR260214", "ERR260174"]
workflow(elasticity_worker, samples, max_procs)
