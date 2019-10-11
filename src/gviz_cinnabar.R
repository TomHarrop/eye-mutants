library(Gviz)
library(GenomicRanges)

options(ucscChromosomeNames=FALSE)

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
cn_cds$cds_name <- "kynurenine 3-monooxygenase"
cn_cds$id <- cn_cds$cds_name

# set up track schemes
set1 <- viridisLite::cividis(9)

set1 <- c("red", "blue", "yellow", "purple")

set1 <- viridis::viridis_pal()(5)

my_scheme <- list(
    shape = "smallArrow",
    background.title = "transparent",
    fontface.title = 1,
    col.title = set1[1],
    col.frame = set1[1],
    thinBoxFeature = c("lincRNA"),
    fontcolor.group = set1[1],
    cex.title = 1,
    cex.group = 0.75
)

grt_scheme <- as.list(c(
    my_scheme,
    col.frame = set1[2],
    fontcolor.group = set1[2],
    fontcolor.title = set1[2],
    size = 1,
    col = set1[2],
    col.line = set1[2],
    fill = set1[2],
    fontface.group = 3,
    transcriptAnnotation = "symbol"
))

cds_scheme <- as.list(c(
    my_scheme,
    col.frame = set1[3],
    fontcolor.group = set1[3],
    col = set1[3],
    col.line = set1[3],
    fill = set1[3],
    fontcolor.title = set1[3],
    fontface.group = 3,
    fill = set1[3]))

aln_scheme <- as.list(c(
    my_scheme,
    fontcolor.title = set1[4]))


cn_chr <- as.character(seqnames(cn_transcripts))[[1]]

cn_start <- min(start(cn_cds))
cn_end <- max(end(cn_cds))

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
displayPars(cn_cdsat) <- cds_scheme




# sequence track
fasta_file <- "data/GCF_003254395.2_Amel_HAv3.1_genomic.fna"
fa <- Biostrings::readDNAStringSet(fasta_file)
names(fa) <- gsub("^([^[:space:]]+).*", "\\1", names(fa)) 
st <- SequenceTrack(fa, chromosome = cn_chr)

# alignment track
at1 <- AlignmentsTrack("output/030_process-aln/Queen_marked.bam",
                       isPaired = TRUE, name = "Queen")
at2 <- AlignmentsTrack("output/030_process-aln/Red1_marked.bam",
                       isPaired = TRUE, name = "Red1")
at3 <- AlignmentsTrack("output/030_process-aln/WT1B3_marked.bam",
                       isPaired = TRUE, name = "WT1B3")
displayPars(at1) <- aln_scheme
displayPars(at2) <- aln_scheme
displayPars(at3) <- aln_scheme


# find the mutated exon
first <- 2760001
last <- 2760027

mut_region <- GRanges(seqnames = cn_chr,
                      ranges = IRanges(start = first, end = last))
exon_hits <- cn_exons[subjectHits(findOverlaps(mut_region, cn_exons))]
smallest_hit <- exon_hits[which.min(width(exon_hits))]

exon_start <- min(start(smallest_hit))
exon_end <- max(end(smallest_hit))

cds_hit <- cn_cds[subjectHits(findOverlaps(mut_region, cn_cds))]
cds_hit$cds_id

grid.newpage()
# full locus
cairo_pdf("test/locus.pdf", width = 8, height = 11)
plotTracks(list(
    cn_grt,
    cn_cdsat,
    at1,
    at2,
    at3,
    st),
           main = "cinnabar homolog locus",
           cex.main = 1,
           fontface.main = 4,
           add = TRUE,
           col.mates = set1[5]
    )
dev.off()


grid.newpage()
# affected exon
cairo_pdf("test/exon.pdf", width = 11, height = 8)
plotTracks(list(cn_cdsat, at1, at2, at3, st),
           from = exon_start - 50,
           to = exon_end + 50,
           main = "Affected exon",
           cex.main = 1,
           fontface.main = 4,
           add = TRUE,
           col.mates = set1[5],
           type = "pileup")

dev.off()




plotTracks(list(at1, at2, at3, st),
           from = 1,
           to = 1e6,
           main = "Affected exon",
           cex.main = 1,
           fontface.main = 4,
           add = TRUE,
           col.mates = "purple")

