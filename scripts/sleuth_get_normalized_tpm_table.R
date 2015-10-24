library(sleuth)

base_dir <- "."

sample_id <- dir(file.path(base_dir,"data/kallisto_quant"))

kal_dirs <- sapply(
    sample_id, 
    function(id) file.path(base_dir, "data/kallisto_quant", id)
)


s2c <- read.table(
    file.path(
        base_dir,
        "data/sleuth/experimental_design.txt"
    ),
    header = TRUE,
    stringsAsFactors=FALSE
) %>%
dplyr::select(
    .,
    sample = run_accession,
    condition
)

so <- sleuth_prep(kal_dirs, s2c, ~ condition)
write.table(
    x = so$obs_norm,
    file = "data/sleuth/tpm_normalized.tsv",
    sep = "\t",
    row.names= FALSE,
    col.names= TRUE,
    quote= FALSE
    )

save(so,
    file= "data/sleuth/sleuth_object.RData")
