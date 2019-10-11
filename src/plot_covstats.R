#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

# dev
covstats_files <- list.files("output/025_pileup/covstats",
                             pattern = ".txt",
                             recursive = TRUE,
                             full.names = TRUE)

hist_files <- list.files("output/025_pileup/hist",
                         pattern = ".txt",
                         recursive = TRUE,
                         full.names = TRUE)


bincov_files <- list.files("output/025_pileup/bincov/",
           pattern = ".txt",
           recursive = TRUE,
           full.names = TRUE)


FreadFiles <- function(file_list){
    names(file_list) <- sub(".txt", "", basename(file_list))
    dt_list <- lapply(file_list, fread)
    dt <- rbindlist(dt_list, idcol = "indiv")
    return(dt)
}

# plot covstats
covstats <- FreadFiles(covstats_files)
ggplot(covstats, aes(x = Ref_GC, y = Avg_fold)) +
    facet_wrap(~indiv) +
    scale_y_log10() +
    geom_point()

# plot hists
hist_data <- FreadFiles(hist_files)
hist_data[]
ggplot(hist_data, aes(x = `#Coverage`, y = numBases)) +
    xlim(c(0, 50)) +
    facet_wrap(~indiv) +
    geom_col(position = "identity")


# plot bincov
# probably look one-per-chr
names(bincov_files) <- sub(".txt", "", basename(bincov_files))
bincov_list <- lapply(bincov_files, fread, skip = 2)
bincov <- rbindlist(bincov_list, idcol = "indiv")
ggplot(bincov, aes(x = RunningPos, y = Cov)) +
    scale_y_log10() +
    facet_wrap(~ indiv) +
    geom_point()


