#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(data.table)
library(GenomicFeatures)
library(VariantAnnotation)

###########
# GLOBALS #
###########


gff_file <- snakemake@input[["gff"]]
vcf_file <- snakemake@input[["vcf"]]
tbi_file <- snakemake@input[["tbi"]]
fa_file <- snakemake@input[["fa"]]
goi_file <- snakemake@input[["goi"]]

coding_file <- snakemake@output[["coding"]]
csd_file <- snakemake@output[["csd"]]
goi_out <- snakemake@output[["goi"]]


# dev
# gff_file <- "data/GCF_003254395.2_Amel_HAv3.1_genomic.gff"
# vcf_file <- "output/050_variant-annotation/variants_filtered.vcf.gz"
# tbi_file <- "output/050_variant-annotation/variants_filtered.vcf.gz.tbi"
# fa_file <- "data/GCF_003254395.2_Amel_HAv3.1_genomic.fna"
# goi_file <- "output/050_variant-annotation/dros_genes.csv"
# coding_file <- "output/050_variant-annotation/coding.Rds"
# csd_file <-"output/050_variant-annotation/csd.vcf"
# goi_out <- "output/050_variant-annotation/goi.vcf"


########
# MAIN #
########

# load the genes of interest
goi <- fread(goi_file,
             fill = TRUE,
             select = 1:4,
             header = TRUE)
goi_ids <- goi[, unique(loc)]

# read the txdb
txdb <- makeTxDbFromGFF(file = gff_file,
                format = "gff3",
                dataSource = "Amel_HAv3.1",
                organism = "Apis mellifera",
                taxonomyId = 7460)

# load the fasta
fa <- FaFile(fa_file, paste0(fa_file, ".fai"))

# read the vcf
tabix_file <- TabixFile(vcf_file, tbi_file)
vcf <- readVcf(tabix_file)

# find interesting changes
# rd <- rowRanges(vcf)
# promoters <- locateVariants(rd, txdb, PromoterVariants())
coding <- predictCoding(vcf, txdb, fa)

# clean up / release mem
rm(vcf, fa, txdb)
gc(TRUE)

# subset by gene-level coding variants
rng <- coding[coding$GENEID %in% goi_ids &
                  coding$CONSEQUENCE != "synonymous"]
vcf_rng <- readVcf(tabix_file, param = rng)

csd <- coding[coding$GENEID == "Csd" & coding$CONSEQUENCE != "synonymous"]
csd_rng <- readVcf(tabix_file, param = csd)

# write the coding changes to file?
saveRDS(coding, coding_file)
writeVcf(vcf_rng, goi_out)
writeVcf(csd_rng, csd_file)

# log
sessionInfo()
