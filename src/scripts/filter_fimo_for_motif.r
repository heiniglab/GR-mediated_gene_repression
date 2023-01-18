suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))

option_list <- list(
  make_option(c("--fimo"),
              type="character",
              help="Path to fimo rds file that should be subsetted"),
  make_option(c("--summit_granges"),
              type="character",
              help="Path to summit granges file that we use to narrow down fimo hits on summitregions"),
  make_option(c("--motif_altname"),
              type="character",
              help="Altname of motif that we filter the fimo results for (case sensitive!)"),
  make_option(c("--outfile"),
              type="character",
              help="Path to subsetted output bed file")
)
opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(library(here, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(dplyr, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(rtracklayer, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(plyranges, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(GenomicRanges, warn.conflicts=F, quietly=T))

fimo_results <- readRDS ( opt$fimo )
ChIPseq_summit_Granges <- readRDS( opt$summit_granges )

# Take 100bp windows around ChIP-seq summits
summit_flank_100bp <- ChIPseq_summit_Granges %>% 
  plyranges::anchor_center() %>% 
  plyranges::mutate(width = 100) 

# narrow the hits of the motif of interest down to the immediate summit region
motifhits <- fimo_results %>% 
  filter(motif_alt_id==opt$motif_altname) %>% 
  filter_by_overlaps(summit_flank_100bp) 

print(paste0("We found " , length(motifhits), " motifhits for ", opt$motif_altname))

export.bed(motifhits, opt$outfile )

