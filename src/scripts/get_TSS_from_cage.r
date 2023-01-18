suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))

option_list <- list(
  make_option(c("--ctss_pool1"), 
              type="character",
              help="Path for ctss file from FANTOM5 of pool 1"),
  make_option(c("--ctss_pool2"), 
              type="character",
              help="Path for ctss file from FANTOM5 of pool 2"),
  make_option(c("--liftoverchain"), 
              type="character",
              help="Path to liftover chain file for mm9 to mm10"),
  make_option(c("--gencode_mm9_geneanno"), 
              type="character",
              help="Path to genomic reference file for mm9"),
  make_option(c("--gencode_mm10_geneanno"), 
              type="character",
              help="GENCODE genomic reference for assembly mm10, prefiltered for gene entries"),
  make_option(c("--outdir"), 
              type="character",
              help="Output directory")
)
opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(library(dplyr, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(BSgenome.Mmusculus.UCSC.mm9, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(ChIPseeker, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(CAGEr, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(rtracklayer, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(GenomicRanges, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(here, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(org.Mm.eg.db, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(TxDb.Mmusculus.UCSC.mm10.knownGene, warn.conflicts=F, quietly=T))

outdir <- here(opt$outdir)

# Workflow is based on the CAGEr vignette
# https://www.bioconductor.org/packages/release/bioc/vignettes/CAGEr/inst/doc/CAGEexp.html


#----------------------   Import CAGE samples from BMDMs
#------------------------------------------------------------------------------------

# get URL for public samples, download and read necessariy columns into ctss file (see snakemake workflow)
# We want to import bone-marrow derived macrophage samples through CAGEr.
# After looking at list of available datasets we can decide what samples best fit our needs.
# Let's see what samples are available through FANTOM5
#data(FANTOM5mouseSamples)
#head(FANTOM5mouseSamples)
# The FANTOM5 dataframe holds descriptions of the samples and the url where they can be retrieved.
# There's an easy way to import samples that match a certain term into a CAGEset object.
#mac_samples <- FANTOM5mouseSamples[grep("macrophage, bone marrow derived",
#                                        FANTOM5mouseSamples[,"description"]),]

print("NOTE: reference genome for the public CAGE files is mm9!")

ce <- CAGEr::CAGEexp( genomeName = "BSgenome.Mmusculus.UCSC.mm9",
                      inputFiles = c(opt$ctss_pool1 ,opt$ctss_pool2),
                      inputFilesType = "ctss",
                      sampleLabels   = c("pool1","pool2")
)
# To actually read in the data into the object we use getCTSS() function, that will add an experiment called tagCountMatrix to the CAGEexp object.
ce <- CAGEr::getCTSS(ce)
ce

#------------------------------------------------------------------------------------
#----------------------   QC
#------------------------------------------------------------------------------------

ncbim37_anno <- rtracklayer::import.gff(opt$gencode_mm9_geneanno)

ce <- annotateCTSS(ce, ncbim37_anno)
colData(ce)[,c("librarySizes", "promoter", "exon", "intron", "unknown")]
plotAnnot(ce, "counts")

corr.m <- plotCorrelation2( ce, samples = "all"
                            , tagCountThreshold = 1, applyThresholdBoth = FALSE
                            , method = "pearson")

#------------------------------------------------------------------------------------
#----------------------   Get read clusters
#------------------------------------------------------------------------------------

print("Merging samples")
#Now we can merge them
ce <- mergeSamples(ce, mergeIndex = c(1,1), 
                   mergedSampleLabels = c("BMDM"))

# redo annotation since this gets reset during merging
ce <- annotateCTSS(ce, ncbim37_anno)

print("The total library size is:")
print(librarySizes(ce))

# Check if data follows a power law distribution
plotReverseCumulatives(ce, fitInRange = c(5, 3000), onePlot = TRUE)

print("Normalizing reads")
# Since we don't really care about making comparisons between different population we could prob just skip the normalization
# The fit range is chosen from the plot. We take the alpha from the ref distribution and set T to a million to get the tag count per million (TPM)
ce <- normalizeTagCount(ce, method = "powerLaw", fitInRange = c(5, 3000), alpha = 1.15, T = 1*10^6)
#mac_CAGEset@tagCountMatrix

print("Cluster the tags")
# After normalization we can cluster the tags.
# Clustering, only seems to work with the CAGEset object (due to some problems with the IRanges column)
# From the CAGEr vignette:
# "Transcription start sites are found in the promoter region of a gene and reflect the transcriptional activity of that promoter (Figure 5). TSSs in the close proximity of each other give rise to a functionally equivalent set of transcripts and are likely regulated by the same promoter elements. Thus, TSSs can be spatially clustered into larger transcriptional units, called tag clusters (TCs) that correspond to individual promoters. CAGEr supports three methods for spatial clustering of TSSs along the genome, two ab initio methods driven by the data itself, as well as assigning TSSs to predefined genomic regions:"

ce <- clusterCTSS(ce,
           threshold=1,
           thresholdIsTpm = TRUE,
           nrPassThreshold = 1,
           method="distclu",
           maxDist=20,
           removeSingletons = TRUE,
           keepSingletonsAbove = 3)

# Let's have a look what the result looks like
head(tagClustersGR(ce, sample = "BMDM"))

# calculate cumulative distribution for every tag cluster in each of the samples
ce <- cumulativeCTSSdistribution(ce, clusters = "tagClusters", useMulticore = T)
# determine the positions of selected quantiles
ce <- quantilePositions(ce, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)

# How many tagclusters do we have in total?
length(tagClustersGR(ce, sample = "BMDM"))

# histogram of interquantile width
plotInterquantileWidth(ce, clusters = "tagClusters", tpmThreshold = 3, qLow = 0.1, qUp = 0.9)

print("Retrieving clusters as GenomicRanges")
clusters_gr <- tagClustersGR(ce, sample="BMDM")

#------------------------------------------------------------------------------------
#----------------------   Liftover coordinates to mm10
#------------------------------------------------------------------------------------

print("Liftover to mm10")

# * Now we can lift over the intervals to mm10
# * Annotate them with peakanno
# * pick the most highly expressed one for each gene

liftover <- function(peaks_gr_mm9){ #input is a GenomicRanges object in mm9 coordinates
  #lift peak locations from mm9 to mm10
  chain <- rtracklayer::import.chain(opt$liftoverchain)
  on.exit( close( file(opt$liftoverchain)) )
  
  peaks_gr_mm10 <- rtracklayer::liftOver(peaks_gr_mm9, chain)
  peaks_gr_mm10 <- GenomicRanges::GRanges(unlist(peaks_gr_mm10))
  
  return(peaks_gr_mm10)
}

mac_cage_mm10 <- liftover( clusters_gr )
ggplot(as.data.frame(mac_cage_mm10), aes(x=width)) +
  geom_histogram(bins = 100)
  
# Liftover coordinates of dominant_ctss
dominant_ctss <- liftover( 
  GRanges(
    seqnames = seqnames(clusters_gr), 
    ranges = IRanges(start = clusters_gr$dominant_ctss,
                     end = clusters_gr$dominant_ctss),
    score=clusters_gr$score)
)

#------------------------------------------------------------------------------------
#-----------------   Annotate TSS clusters to reference gene coordinates
#------------------------------------------------------------------------------------
print("Annotate TSS coordinates")
# Use coordinates of the dominant ctss downstream

mac_cage_anno <- ChIPseeker::annotatePeak(dominant_ctss, 
                                          tssRegion=c(-1000, 1000), #more stringent than default
                                          level = "gene",
                                          TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, 
                                          annoDb = "org.Mm.eg.db")

# For those that are reasonably close to a TSS,
# check for each gene, which position has the highest score.

mac_cage_maxscore <- as.data.frame(mac_cage_anno) %>% 
  filter(abs(distanceToTSS)<=30000) %>%
  mutate(SYMBOL=as.factor(SYMBOL)) %>%
  filter(SYMBOL!="") %>%
  group_by(SYMBOL)%>%
  filter(score == max(score))%>%
  filter(distanceToTSS == min(distanceToTSS )) # for tied score, use shorter distance

nrow(mac_cage_maxscore)

ggplot(mac_cage_maxscore, aes(x=distanceToTSS)) +
  geom_histogram(bins = 100)

#------------------------------------------------------------------------------------
#----------   retrieve gene coordinates and promoterregion from reference
#------------------------------------------------------------------------------------

gencode_mm10_geneanno <- rtracklayer::import.gff(opt$gencode_mm10_geneanno)
genecoords <- as.data.frame(gencode_mm10_geneanno) %>%
  dplyr::select("seqnames","start","end","strand","gene_id") %>%
  mutate(score=0) %>%
  dplyr::mutate(gene_id=gsub("\\.[0-9]*$","",gene_id)) %>%
  dplyr::filter(!"gene_id"=="")%>%
  dplyr::select("seqnames","start","end","gene_id","score","strand")


gencode_mm10_promoterregions <- promoters(gencode_mm10_geneanno)
gencode_mm10_promoterregions <- as.data.frame(gencode_mm10_promoterregions) %>%
  dplyr::select("seqnames","start","end","strand","gene_id") %>%
  mutate(score=0) %>%
  dplyr::mutate(gene_id=gsub("\\.[0-9]*$","",gene_id)) %>%
  dplyr::filter(!"gene_id"=="")%>%
  dplyr::select("seqnames","start","end","gene_id","score","strand")
  
#------------------------------------------------------------------------------------
#-----------------   export files
#------------------------------------------------------------------------------------

write.table(genecoords, 
            file = paste0(outdir,"reference_genecoords.bed"),
            sep="\t", 
            col.names = FALSE,
            quote=FALSE,
            row.names = FALSE)

write.table(gencode_mm10_promoterregions, 
            file = paste0(outdir,"reference_promoterregions.bed"),
            sep="\t", 
            col.names = FALSE,
            quote=FALSE,
            row.names = FALSE)


rtracklayer::export.bed(as.data.frame(mac_cage_mm10),
                        paste0(outdir, "mac_cage_tssclusterregions.bed"),
                        format="bed")

rtracklayer::export.bed(as.data.frame(dominant_ctss),
                        paste0(outdir, "mac_cage_dominant_ctss.bed"),
                        format="bed")

rtracklayer::export.bed(mac_cage_maxscore,
                        paste0(outdir,"mac_cage_maxscore.bed"),
                        format="bed")
