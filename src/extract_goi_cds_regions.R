#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output", append = TRUE)

library(data.table)
library(GenomicFeatures)

gff_file <- snakemake@input[["gff"]]
coding_file <- snakemake@input[["coding"]]
assoc_file <- snakemake@input[["assoc"]]
regions_file <- snakemake@output[["regions"]]

# dev
# gff_file <- "data/GCF_003254395.2_Amel_HAv3.1_genomic.gff"
# coding_file <- "output/050_variant-annotation/coding.Rds"
# assoc_file <- "output/060_plink/goi.assoc"

# make a TXDB
txdb <- makeTxDbFromGFF(file = gff_file,
                        format = "gff3",
                        dataSource = "Amel_HAv3.1",
                        organism = "Apis mellifera",
                        taxonomyId = 7460)
# get features
all_cds <- cdsBy(txdb, "gene")
all_tx <- transcriptsBy(txdb, "gene")
all_exons <- exonsBy(txdb, "gene")

# which SNP do we want to look at
assoc <- fread(assoc_file)
setorder(assoc, P, na.last = TRUE)
snpoi <- assoc[P == min(P, na.rm = TRUE), c(SNP)]

# which gene is that SNP in
coding <- readRDS(coding_file)
goi <- unique(coding[snpoi]$GENEID)

# find the transcripts that contain the snp
goi_tx <- all_tx[[goi]]
tx_hits <- subjectHits(
    findOverlaps(coding[snpoi], goi_tx, type = "within"))
overlapping_tx <- unique(goi_tx[tx_hits])

# find the longest transcript
overlapping_tx[which.max(width(overlapping_tx))]$tx_name

# find the longest overlapping cds
goi_cds <- all_cds[[goi]]
cds_hits <- subjectHits(
    findOverlaps(coding[snpoi], goi_cds, type = "within"))
overlapping_cds <- unique(goi_cds[cds_hits])
longest_cds <- as.data.table(overlapping_cds)[
    , sum(width), by = cds_name][which.max(V1), cds_name]

# extract the data for the longest cds exons
cdsoi <- goi_cds[goi_cds$cds_name == longest_cds]
cds_dt <- as.data.table(cdsoi)

# generate the samtools region format
outlines <- cds_dt[, paste0(seqnames, ":", start, "-", end)]
writeLines(outlines, regions_file, sep = " ")

# log
sessionInfo()


