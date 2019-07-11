"""Calculate the elasticities for a set of built models."""

from os.path import isfile
from micom import load_pickle
from micom.elasticity import exchange_elasticities
from micom.workflows import workflow


max_procs = 20


def elasticities(sam):
    """Get the exchange elasticities for a sample."""
    model_file = "models/" + sam + ".pickle"
    out_file = "data/elasticities_" + sam + ".csv"
    if isfile(out_file):
        return
    com = load_pickle(model_file)
    elast = exchange_elasticities(com, fraction=0.5, progress=False)
    elast.to_csv(out_file, index=False)


samples = ["ERR260275", "ERR260214", "ERR260174"]
workflow(elasticities, samples, max_procs)
