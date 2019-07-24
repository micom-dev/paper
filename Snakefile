samples = ["ERR260275", "ERR260214", "ERR260174"]

rule all:
    input:
        "data/growth_rates.csv",
        "data/knockouts.csv",
        "figures/gcs.svg",
        "figures/circos.svg",
        "figures/media.png",
        expand("figures/elasticities_{s}.png", s=samples)

rule collapse:
    input:
        "data/abundances.csv"
    output:
        "data/species.csv"
    threads: 1
    script:
        "workflows/species.py"

rule build_models:
    input:
        "data/species.csv",
        "data/western_diet.csv"
    output:
        directory("data/models")
    threads: 32
    script:
        "workflows/build_models.py"

rule tradeoff:
    input:
        "data/models",
        "data/recent.csv"
    output:
        "data/tradeoff.csv"
    threads: 32
    script:
        "workflows/tradeoff.py"

rule media_and_rates:
    input:
        "data/models"
    output:
        "data/growth_rates.csv",
        "data/minimal_imports.csv",
        "data/minimal_fluxes.csv.gz"
    threads: 32
    script:
        "workflows/media_and_gcs.py"

rule knockouts:
    input:
        "data/models"
    output:
        "data/knockouts.csv"
    threads: 32
    script:
        "workflow/knockouts.py"

rule elasticities:
    input:
        "data/models"
    output:
        expand("data/elasticities_{s}.csv", s=samples)
    threads: 3
    script:
        "workflow/elasticities.py"

rule rate_figures:
    input:
        "data/replication_rates.csv",
        "data/tradeoff.csv"
    output:
        "figures/dists.svg",
        "figures/non_zero.svg",
        "figures/percent_growing.svg",
        "figures/gcs.svg",
        "figures/community_growth.svg",
        "figures/rate_vs_abundance.png"
    threads: 1
    script:
        "workflow/growth_rate_figs.py"

rule knockout_figures:
    input:
        "data/knockouts.csv"
    output:
        "figures/circos.svg"
    threads: 1
    script:
        "workflows/knockout_figs.py"

rule exchange_figures:
    input:
        "data/growth_rates.csv",
        "data/minimal_imports.csv",
        "data/minimal_fluxes.csv.gz",
        "data/metabolites.csv"
    output:
        "figures/media.png",
        "figures/scfas_prod.svg",
        "figures/scfas_consumption.svg",
        "figures/scfas_net.svg",
        "figures/scfas.svg",
        "figures/individual_media.png"
    threads: 1
    script:
        "workflows/exchange_figs.py"

rule elasticity_figures:
    input:
        expand("data/elasticities_{s}.csv", s=samples)
    output:
        expand("figures/elasticities_{s}.png", s=samples)
    threads: 1
    script:
        "workflows/elasticity_figs.py"
