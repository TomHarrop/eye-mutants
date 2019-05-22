#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)

key_file <- snakemake@input[["key"]]
pheno_file <- snakemake@output[["pheno"]]

# read data
key_data <- fread(key_file)

# make the plink table
plink_pheno <- copy(key_data[caste == "drone", .(
    Family = our_id,
    Individual = our_id,
    Phenotype = ifelse(grepl("^WT", our_id), 0, 1))])

# match plink's weird expectations on underscores
plink_pheno[grep("\\s", Family),
            c("Family", "Individual") := tstrsplit(Family, "\\s")]

# write output
fwrite(plink_pheno, pheno_file, col.names = FALSE, sep = " ")

# log session
sessionInfo()
