suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))

option_list <- list(
  make_option(c("--chipseq_bam"),
              type="character",
              help="Path to BAM file with deduplicated reads from DexLPS condition"),
  make_option(c("--chipseq_summits"),
              type="character",
              help="Path to summit file of IDR peaks"),
  make_option(c("--nr3c1fullsitematches"),
              type="character",
              help="Path to homer hits of nr3c1 fullsite motif"),
  make_option(c("--nr3c1halfsitematches"),
              type="character",
              help="Path to homer hits of nr3c1 halfsite motif"),
  make_option(c("-o", "--outdir"),
              type="character",
              help="Path to output directory"))

opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(library(dplyr, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(rtracklayer, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(plyranges, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(genomation, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(here, warn.conflicts=F, quietly=T))

#------ chipseq summits
#----------------------------------------------
chipseq_summits <- read.table(opt$chipseq_summits)
chipseq_summits <- GRanges(seqnames = chipseq_summits[,c("V1")],
                           ranges = IRanges(start=chipseq_summits[,c("V2")],
                                            end=chipseq_summits[,c("V3")]-1)) # to make up for 0 vs 1 encoding
chipseq_summits$id <- c(1:length(chipseq_summits))

chipseq_summitranges <- chipseq_summits %>% 
  plyranges::anchor_center() %>% 
  plyranges::mutate(width = 1000) 

#------ motifmatches
#----------------------------------------------

nr3c1fullsitematches <- rtracklayer::import.bed(opt$nr3c1fullsitematches)
# remove weird repetition of seqname
seqlevels(nr3c1fullsitematches) <- gsub(" .*$","",seqlevels(nr3c1fullsitematches))
nr3c1fullsitematches <- GenomeInfoDb::keepStandardChromosomes(nr3c1fullsitematches, pruning.mode = "tidy")
names(nr3c1fullsitematches) <- c(1:length(nr3c1fullsitematches))
nr3c1fullsitematches_ranges <- nr3c1fullsitematches %>% 
  plyranges::anchor_center() %>% 
  plyranges::mutate(width = 1000) 

#------ make subsets by combining chipseq summits and motifhits
#----------------------------------------------
# only use the4 summitranges containing an nr3c1 motif
chipseq_summitranges_inner_hits <- chipseq_summitranges %>% 
  plyranges::join_overlap_inner(nr3c1fullsitematches)
width(chipseq_summitranges_inner_hits)
length(chipseq_summitranges_inner_hits)

# only use the nr3c1 coordinates that fall within summitranges
chipseq_summitranges_intersect_hits <- chipseq_summitranges %>% 
  plyranges::join_overlap_intersect(nr3c1fullsitematches)
width(chipseq_summitranges_intersect_hits)
length(chipseq_summitranges_intersect_hits)
# size it up from 14bp to 1000bp around the motifhit
chipseq_summitranges_intersect_hits <- 
  chipseq_summitranges_intersect_hits %>%
  plyranges::anchor_center() %>% 
  plyranges::mutate(width = 1000) 


#----------------------------------------------
#------ compute score matrix
#----------------------------------------------

# for each summitregion, get coverage of how many reads in bam file overlap it at each position
# aggregate it across all peaks

sm_summitranges <- ScoreMatrix(target = opt$chipseq_bam,
                               windows = chipseq_summitranges,
                               weight.col = "score")

sm_nr3c1fullsitematches_ranges <- ScoreMatrix(target = opt$chipseq_bam,
                                              windows = nr3c1fullsitematches_ranges,
                                              weight.col = "score")

sm_summitranges_w_nr3c1fullsitehit <- ScoreMatrix(target = opt$chipseq_bam,
                                                  windows = chipseq_summitranges_inner_hits,
                                                  weight.col = "score")

sm_nr3c1fullsitehits_within_summitranges <- ScoreMatrix(target = opt$chipseq_bam,
                                                        windows = chipseq_summitranges_intersect_hits,
                                                        weight.col = "score")

test_coverage <- chipseq_summitranges_intersect_hits %>% plyranges::compute_coverage()
score(test_coverage)
#----------------------------------------------
#------ export data
#----------------------------------------------

saveRDS(sm_summitranges,
        paste0(opt$outdir,"sm_summitranges.rds"))
saveRDS(sm_nr3c1fullsitematches_ranges,
        paste0(opt$outdir,"sm_nr3c1fullsitematches_ranges.rds"))
saveRDS(sm_summitranges_w_nr3c1fullsitehit,
        paste0(opt$outdir,"sm_summitranges_w_nr3c1fullsitehit.rds"))
saveRDS(sm_nr3c1fullsitehits_within_summitranges,
        paste0(opt$outdir,"sm_nr3c1fullsitehits_within_summitranges.rds"))
