suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))

option_list <- list(
  make_option(c( "--ABC_DexLPS_all"),
              type="character",
              help="Path to file with all ABC predictions passing 0.02 (in DexLPS condition)"),
  make_option(c("--ABC_LPS_all"),
              type="character",
              help="Path to file with all ABC predictions passing 0.02 (in LPS condition)"),
  make_option(c("--fimo_results_summitregion"),
              type="character",
              help="Path to rds file of fimo motifcounts within summitregions"),
  make_option(c("--fimo_results_dexlps"),
              type="character",
              help="Path to rds file of fimo motifcounts within ABC regions (in DexLPS condition)"),
  make_option(c("--fimo_results_lps"),
              type="character",
              help="Path to rds file of fimo motifcounts within ABC regions (in LPS condition)"),
  make_option(c( "--chipseq_ranges"),
              type="character",
              help="Path to summit file of IDR peaks"),
  make_option(c("--outdir"),
              type="character",
              help="Path to output directory")
)

opt <- parse_args(OptionParser(option_list=option_list))

dir.create(opt$outdir)

#change default for stringAsFactors
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages(library(dplyr, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(here, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(plyranges, warn.conflicts=F, quietly=T))

#set defaults for ggplot2 figures
theme_update(panel.background = element_rect(fill = "transparent", colour = NA),
             plot.background = element_rect(fill = "transparent", colour = NA),
             legend.background = element_rect(fill = "transparent", colour = NA),
             legend.key = element_rect(fill = "transparent", colour = NA),
             text=element_text(size=10, family = "ArialMT", colour="black"),
             title=element_text(size=10, family="ArialMT", colour="black"),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             axis.text = element_text(size=10, family="ArialMT", colour="black"),
             axis.line = element_line(colour="black"),
             axis.ticks = element_line(colour="black"))

set.seed(12345)

#-------------------------------
## read in  data
#-------------------------------
ABC_DexLPS_all <- read.delim(opt$ABC_DexLPS_all) %>% plyranges::as_granges(., seqnames=chr)
ABC_LPS_all <- read.delim(opt$ABC_LPS_all) %>% plyranges::as_granges(., seqnames=chr)
fimo_results_dexlps <- readRDS(opt$fimo_results_dexlps)
fimo_results_lps <- readRDS(opt$fimo_results_lps)
fimo_results_summitregion <- readRDS(opt$fimo_results_summitregion)
chipseq_ranges <- readRDS(opt$chipseq_ranges)

# assign a name to allow for matching for prox based assignment later
chipseq_ranges$name <- paste(seqnames(chipseq_ranges),
                             start(chipseq_ranges),
                             end(chipseq_ranges),
                             sep="_")

#--------------------------------------------
#--------------------------------------------
## prepare motifcounts
#--------------------------------------------
#--------------------------------------------

## summitregions
#--------------------------------------------

get_summitregion_motifcounts <- function(summits, fimo_results){
  
  # Take 100bp windows around ChIP-seq summits to recreate the original query
  fimo_queries_summitregion <- summits %>% 
    plyranges::anchor_center() %>% 
    plyranges::mutate(width = 100)
  
  ## intersect queries with fimo hits
  summitregion_leftjoin_query_hits <- fimo_queries_summitregion %>% 
    plyranges::join_overlap_left(fimo_results)
  
  # aggregate motifcounts per query
  # seqnames is needed for the train test split later
  summitregion_leftjoin_query_hits_motifsaggregated <-
    as.data.frame(summitregion_leftjoin_query_hits) %>% 
    group_by(name, seqnames, motif_alt_id) %>%
    summarize(motifcount = n())
  # if a region has no motifmatches at all, it get's an NA which showes up as motifname after doing pivor_wider
  
  # cast it in a way, so we have unique regions as rows and all observed motifs as columns
  motifcounts <- 
    summitregion_leftjoin_query_hits_motifsaggregated %>%
    tidyr::pivot_wider(names_from = motif_alt_id,
                       values_from = motifcount,
                       values_fill = 0) %>%
    dplyr::select(!'NA') # remove NA motif that got introduced by region without any matches
  
  return(motifcounts)
}

motifcounts_summitregion <- get_summitregion_motifcounts(
  chipseq_ranges,
  fimo_results_summitregion
)

## ABC enhancerregions
#--------------------------------------------

get_ABC_motifcounts <- function(ABC_results, fimo_results){
  # intersect queries with fimo hits
  ABC_results_unique <- ABC_results %>% unique()
  ABC_results_leftjoin_query_hits <- ABC_results_unique %>% 
    plyranges::join_overlap_left(fimo_results)
  
  # aggregate motifcounts per query
  ABC_results_leftjoin_query_hits_motifsaggregated <-
    as.data.frame(ABC_results_leftjoin_query_hits) %>% 
    group_by(name, seqnames ,motif_alt_id) %>%
    summarize(motifcount = n())
  
  # cast it in a way, so we have unique regions as rows and all observed motifs as columns
  motifcounts_abcregion <- 
    ABC_results_leftjoin_query_hits_motifsaggregated %>%
    tidyr::pivot_wider(names_from = motif_alt_id,
                       values_from = motifcount,
                       values_fill = 0) %>%
    dplyr::select(!'NA') # remove NA motif that got introduced by region without any matches
  
  return(motifcounts_abcregion)
}

motifcounts_abcregion_dexlps <- get_ABC_motifcounts(
  ABC_DexLPS_all,
  fimo_results_dexlps
)
motifcounts_abcregion_lps <- get_ABC_motifcounts(
  ABC_LPS_all,
  fimo_results_lps
)

#--------------------------------------------
## ASSIGNMENTS
#--------------------------------------------

## prox based
#--------------------------------------------

# these are mgi_symbols, the ones from ABC are ensemble IDs. Is this a problem?
assignment_summit_prox <- as.data.frame(chipseq_ranges) %>%
  dplyr::rename(anno=mgi_symbol) %>%
  dplyr::select(name, anno)

## hybrid
#--------------------------------------------
# we need an assignment of the regionIDs (from the ChIPseq summtiregion)
# to the ABC derived assignments

get_assignment_summit_ABCregion <- function(summit, ABCregions){
  # what summit lies in what ABC region
  assignment_summit_ABCregion <- summit %>% 
    plyranges::join_overlap_left(ABCregions)
  
  assignment_summit_ABCregion <- as.data.frame(assignment_summit_ABCregion) %>% 
    dplyr::select(name.x, TargetGene, ABC.Score, ABC.Score.Numerator, class, isSelfPromoter) %>%
    magrittr::set_colnames(c("name","anno","abcscore","abcnumerator","class","isSelfPromoter")) 
  
  return(assignment_summit_ABCregion)
}

assignment_summits_ABCregion_dexlps <- get_assignment_summit_ABCregion(
  chipseq_ranges,
  ABC_DexLPS_all
)
assignment_summits_ABCregion_lps <- get_assignment_summit_ABCregion(
  chipseq_ranges,
  ABC_LPS_all
)

## ABC based
#--------------------------------------------
assignment_abcregion_dexlps <- as.data.frame(ABC_DexLPS_all) %>% 
  dplyr::select(name,TargetGene,ABC.Score,ABC.Score.Numerator,class,isSelfPromoter) %>%
  magrittr::set_colnames(c("name","anno","abcscore","abcnumerator","class","isSelfPromoter")) %>%
  dplyr::mutate(anno=gsub("\\.[0-9]*$","",anno))

assignment_abcregion_lps <- as.data.frame(ABC_LPS_all) %>% 
  dplyr::select(name,TargetGene,ABC.Score,ABC.Score.Numerator,class,isSelfPromoter) %>%
  magrittr::set_colnames(c("name","anno","abcscore","abcnumerator","class","isSelfPromoter")) %>%
  dplyr::mutate(anno=gsub("\\.[0-9]*$","",anno))

#--------------------------------------------
## export objects
#--------------------------------------------

# assignments
saveRDS(assignment_summit_prox, paste0(opt$outdir,"assignment_summit_prox.rds"))
saveRDS(assignment_summits_ABCregion_dexlps, paste0(opt$outdir,"assignment_summits_ABCregion_dexlps.rds"))
saveRDS(assignment_summits_ABCregion_lps, paste0(opt$outdir,"assignment_summits_ABCregion_lps.rds"))
saveRDS(assignment_abcregion_dexlps, paste0(opt$outdir,"assignment_abcregion_dexlps.rds"))
saveRDS(assignment_abcregion_lps, paste0(opt$outdir,"assignment_abcregion_lps.rds"))
# motifcoutns
saveRDS(motifcounts_summitregion, paste0(opt$outdir,"motifcounts_summitregion.rds"))
saveRDS(motifcounts_abcregion_dexlps, paste0(opt$outdir,"motifcounts_abcregion_dexlps.rds"))
saveRDS(motifcounts_abcregion_lps, paste0(opt$outdir,"motifcounts_abcregion_lps.rds"))
