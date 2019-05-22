library(Gviz)

gff_file <- "data/GCF_003254395.2_Amel_HAv3.1_genomic.gff"

txdb_file <- "test/txdb.sqlite"
if(!file.exists(txdb_file)) {
    txdb <- GenomicFeatures::makeTxDbFromGFF(gff_file)
    AnnotationDbi::saveDb(txdb, txdb_file)
} else {
    txdb <- AnnotationDbi::loadDb(txdb_file)
}


# set up transcripts, exons and coding sequences
tx <- transcriptsBy(txdb, "gene")
exons <- exonsBy(txdb, "gene")
cds <- cdsBy(txdb, "gene")

# find cinnabar in the txdb
cn_transcripts <- tx$LOC551858
cn_exons <- exons$LOC551858
cn_cds <- cds$LOC551858
cn_cds$cds_name <- "cinnabar"
cn_cds$id <- cn_cds$cds_name

# set up tracks
set1 <- RColorBrewer::brewer.pal(9, "Set1")
my_scheme <- list(
    shape = "smallArrow",
    background.title = "transparent",
    fontface.title = 1,
    col.title = "black",
    col.frame = "black",
    col = "black",
    col.line = "black",
    thinBoxFeature = c("lincRNA"),
    fontcolor.group = "black",
    cex.title = 1,
    cex.group = 0.75
)

grt_scheme <- as.list(c(
    my_scheme,
    fontcolor.title = set1[1],
    size = 1,
    fill = set1[1],
    fontface.group = 3,
    transcriptAnnotation = "symbol"
))

cds_scheme <- as.list(c(
    my_scheme,
    fontcolor.title = set1[3],
    fontface.group = 3,
    fill = set1[3]))

cn_chr <- as.character(seqnames(cn_transcripts))[[1]]

cn_start <- min(start(cn_cds))
cn_end <- max(end(cn_cds))
options(ucscChromosomeNames=FALSE)
cn_grt <- GeneRegionTrack(txdb,
                          chromosome = cn_chr,
                          start = cn_start,
                          end = cn_end,
                          name = "NCBI transcripts")
displayPars(cn_grt) <- grt_scheme

cn_cdsat <- AnnotationTrack(cn_cds,
                            name = "CDS",
                            group = "cds_name",
                            groupAnnotation = "id")

# sequence track
fasta_file <- "data/GCF_003254395.2_Amel_HAv3.1_genomic.fna"
fa <- readDNAStringSet(fasta_file)
names(fa) <- gsub("^([^[:space:]]+).*", "\\1", names(fa)) 
st <- SequenceTrack(fa, chromosome = cn_chr)

# alignment track
at1 <- AlignmentsTrack("output/030_process-aln/Queen_marked.bam",
                       isPaired = TRUE, name = "Queen")
at2 <- AlignmentsTrack("output/030_process-aln/Red1_marked.bam",
                       isPaired = TRUE, name = "Red1")
at3 <- AlignmentsTrack("output/030_process-aln/WT1B3_marked.bam",
                       isPaired = TRUE, name = "WT1B3")


first <- 2760001
last <- 2760027
middle <-  first + ((last - first) %/% 2)
offset <- 600

grid.newpage()
plotTracks(list(cn_grt, cn_cdsat, at1, at2, at3, st),
           #from = middle - offset,
           #to = middle + offset,
           main = "Cpr",
           cex.main = 1,
           fontface.main = 4,
           add = TRUE,
           col.mates = "purple")


