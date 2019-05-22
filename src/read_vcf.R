library(vcfR)
library(data.table)
library(ggplot2)

vcf_file <- "output/040_freebayes/variants_filtered.vcf"
dna_file <- "data/GCF_003254395.2_Amel_HAv3.1_genomic.fna"
ann_file <- "data/GCF_003254395.2_Amel_HAv3.1_genomic.gff"

vcf <- read.vcfR(vcf_file)
dna <- ape::read.dna(dna_file, format = "fasta")

gff_table_all <- read.table(ann_file,
                            sep = "\t",
                            quote = "",
                            stringsAsFactors = FALSE)

keep_features = c("exon")

gff_filtered <- gff_table_all[gff_table_all$V3 %in% keep_features &
                                  gff_table_all$V1 == "NC_037640.1", ]


chrom <- create.chromR(name = "NC_037640.1",
                       vcf = vcf,
                       seq = dna[grepl("NC_037640.1", names(dna))],
                       ann = gff_filtered)

chrom_masked <- masker(chrom, min_QUAL = 20)
chrom_processed <- proc.chromR(chrom, verbose = TRUE)


plot(chrom)
plot(chrom_masked)
plot(chrom_processed)

chromoqc(chrom_masked, xlim = c(11771679-100, 11781139+100))

x <- data.table(vcf@fix)
x[, QUAL := as.numeric(QUAL)]
x[, POS := as.numeric(POS)]

ggplot(x, aes(x = POS, y = QUAL)) +
    xlim(c(11771679-100, 11781139+100)) +
    scale_y_log10() +
    geom_hline(yintercept = 100) +
    geom_point()


