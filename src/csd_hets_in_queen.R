library(data.table)
library(GenomicFeatures)
library(VariantAnnotation)


coding <- readRDS("output/050_variant-annotation/coding.Rds")
vcf_file <- "output/050_variant-annotation/variants_filtered.vcf.gz"
tbi_file <- "output/050_variant-annotation/variants_filtered.vcf.gz.tbi"

tabix_file <- TabixFile(vcf_file, tbi_file)

# subset by gene-level coding variants
csd <- coding[coding$GENEID == "Csd" & coding$CONSEQUENCE != "synonymous"]
csd_rng <- readVcf(tabix_file, param = csd)

# only look at *coding* variants that are heterozygous in the queen
csd_dt <- data.table(geno(csd_rng)$"GT", keep.rownames = TRUE)
keep_loci <- csd_dt[Queen %in% c("0/1", "1/0"), unique(rn)]

csd_hets <- csd_rng[keep_loci]

writeVcf(csd_hets, "test/csd_hets.vcf")
