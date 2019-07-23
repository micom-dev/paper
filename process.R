library(mbtools)
library(futile.logger)
library(drake)

flog.appender(appender.tee("mbtools.log"))
options(mc.cores = 20)
pattern <- "(\\w+)_(\\d+)\\.fastq.gz"
annotations <- c("id", "direction")


save_tables <- function(sl, rates) {
    fwrite(sl$abundance, "data/abundances.csv")
    fwrite(rates$rate, "data/replication_rates.csv")
}


plan <- drake_plan(
    files = find_read_files(
        "data/raw",
        pattern, annotations
    ),
    qual = quality_control(files, min_score = 20),
    plot_qual = ggsave(
        file_out("figures/quals.png"),
        plot = qual$quality_plot + theme_minimal(),
        width = 10,
        height = 4,
        dpi = 200),
    processed = preprocess(
        files,
        out_dir = "data/preprocessed",
        trimLeft = 10,
        maxEE = 10),
    aligned = align_short_reads(
        processed,
        reference = "../refs/ABVF_SP_genomes.fna.gz",
        alignment_dir = "data/alignments"
    ),
    sl = slimm(
        aligned,
        reports = file_out("data/slimm_reports"),
        database = file_in("../refs/ABVF_SP_CMP_genomes.sldb"),
        bin_width = 100
    ),
    rates = replication_rates(sl),
    tables = save_tables(sl, rates)
)

make(plan, lock_envir = FALSE)
