suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))

option_list <- list(
  make_option(c( "--summit_granges"),
              type="character",
              help="Path to rds file of summits in granges format with directionschange and distancetoTSS as additional metadata columns"),
  make_option(c("--memedb_expressed"),
              type="character",
              help="Path to memedb file filtered for motifs where TFs are expressed in 4sU"))

opt <- parse_args(OptionParser(option_list=option_list))
  
suppressPackageStartupMessages(library(here, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(memes, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(universalmotif, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(dplyr, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(plyranges, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(GenomicRanges, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(BSgenome.Mmusculus.UCSC.mm10, warn.conflicts=F, quietly=T))

#-------------------------------
## Import reference for sequence
#-------------------------------
mm.genome <- BSgenome.Mmusculus.UCSC.mm10

#-------------------------------
## read in prepared data
#-------------------------------

ChIPseq_summit_Granges <- readRDS(opt$summit_granges)

# Take 100bp windows around ChIP-seq summits
summit_flank_100bp <- ChIPseq_summit_Granges %>% 
  plyranges::anchor_center() %>% 
  plyranges::mutate(width = 100) 

# Take 100bp windows around ChIP-seq summits
summit_flank_1000bp <- ChIPseq_summit_Granges %>% 
  plyranges::anchor_center() %>% 
  plyranges::mutate(width = 1000) 

meme_db_expressed <- readRDS(opt$memedb_expressed)


# to_list() converts the database back from data.frame format to a standard `universalmotif` object.
options(meme_db = to_list(meme_db_expressed, extrainfo = FALSE))

# where is meme installed 
my_memepath="~/software/meme/bin/"
check_meme_install(meme_path=my_memepath)

summit_flank_100bp_seq <- summit_flank_100bp %>%
  get_sequence(mm.genome)

summit_flank_1000bp_seq <- summit_flank_1000bp %>%
  get_sequence(mm.genome)


#-------------------------------
## run fimo
#-------------------------------

fimo_results <-
  runFimo(summit_flank_1000bp_seq,
          meme_db_expressed,
          meme_path=my_memepath)

saveRDS (fimo_results, here("results/current/memes_bioc/fimo_1000bp/fimo.rds") )

print("Finished running fimo")
print("Analysis DONE")
