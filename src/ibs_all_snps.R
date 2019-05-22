library(SNPRelate)
library(data.table)

vcf_file <- "output/040_freebayes/variants_filtered.vcf"
gds_file <- tempfile(fileext = ".gds")

# make a gds dataset
snpgdsVCF2GDS(vcf_file, gds_file, method = "copy.num.of.ref" )
gds_data <- snpgdsOpen(gds_file)

# select SNPs (use all for testing)
snp_set <- snpgdsSelectSNP(gds_data,
                           maf = 0.1,
                           autosome.only = FALSE,
                           verbose = TRUE,
                           missing.rate = 0.2
)

snp_set_ids <-  unlist(snp_set)

# Get samples where missing rate is higher than sample missing quantile
sample_missing_rates <- snpgdsSampMissRate(gdsobj = gds_data,
                                           snp.id = snp_set_ids,
                                           with.id = TRUE)

kept_indivs <- sample_missing_rates[
    sample_missing_rates < 0.8]

# run identity by state
ibs_res <- snpgdsIBS(gds_data,
                     sample.id = names(kept_indivs),
                     snp.id = snp_set_ids,
                     autosome.only = FALSE)

# convert to dist for clustering
ibs_df <- data.frame(ibs_res$ibs, row.names = ibs_res$sample.id)
colnames(ibs_df) <- ibs_res$sample.id

# cluster the ibs results
d <- as.dist(1 - ibs_df)
hc <- hclust(d, method = "ward.D2")
indiv_order <- hc$labels[hc$order]

# generic plot of the clusters
xw <- grid::convertUnit(unit(210 - 20, "mm"), "in", valueOnly = TRUE)
cairo_pdf("dendro.pdf", pointsize = 8, width = xw, height = xw)
plot(hc)
dev.off()

# data.table of ibs results for plotting
ibs_dt <- data.table(ibs_df, keep.rownames = TRUE)
setnames(ibs_dt, "rn", "i1")
ibs_pd <- melt(ibs_dt, id.vars = "i1", variable.name = "i2")

ibs_pd[, i1 := factor(i1, levels = indiv_order)]
ibs_pd[, i2 := factor(i2, levels = indiv_order)]

# set up label colours
label_types <- gsub("[[:digit:]]+", "", indiv_order)
label_cols <- plyr::mapvalues(label_types,
                              from = unique(label_types),
                              to = RColorBrewer::brewer.pal(length(unique(label_types)), "Set1"))

# plot heatmap of similarity
gp <- ggplot(ibs_pd, aes(x = i1, y = i2, fill = value)) +
    theme_grey(base_size = 8) +
    theme(axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.5,
                                     colour = label_cols),
          axis.text.y = element_text(colour = label_cols)) +
    xlab(NULL) + ylab(NULL) + coord_fixed() +
    scale_fill_viridis_c(guide = guide_colorbar(title = "IBS"),
                         limits = c(0,1)) +
    geom_raster()

ggsave("test.pdf",gp, width = 210 - 20, height = 210 - 20, units = "mm")

