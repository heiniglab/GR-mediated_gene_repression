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
suppressPackageStartupMessages(library(memes, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(universalmotif, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(GenomicRanges, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(BSgenome.Mmusculus.UCSC.mm10, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(plyranges, warn.conflicts=F, quietly=T))

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

#-------------------------------
# define inputs
#-------------------------------

summit_flank_100bp_seq <- summit_flank_100bp %>%
  get_sequence(mm.genome)

summit_flank_1000bp_seq <- summit_flank_1000bp %>%
  get_sequence(mm.genome)

summit_flank_seq_bydirchange <- summit_flank_100bp %>%
  # remove unchanged ones and only compare "up" vs "down"
  filter(directionchange !="ns")%>%
  # Get a list of chip peaks belonging to each set
  split(mcols(.)$directionchange) %>%
  # look up the DNA sequence of each peak within each group
  get_sequence(mm.genome)

#-------------------------------
## up vs downregulation
# run by directionchange to discover consensus motif separately
#-------------------------------

#-------------------------------
# STREME
#-------------------------------

print("Start running streme for 100bp")

stremeout_100bp_down <-  here("results/current/memes_bioc/streme_100bp_down/streme.xml")
if (!file.exists( stremeout_100bp_down )){
  runStreme(summit_flank_seq_bydirchange[["down"]], control="shuffle", objfun="de",
            meme_path="~/software/meme/bin/", silent=FALSE,
            outdir = dirname(stremeout_100bp_down))
}

stremeout_100bp_up <- here("results/current/memes_bioc/streme_100bp_up/streme.xml")
if (!file.exists( stremeout_100bp_up )){
  runStreme(summit_flank_seq_bydirchange[["up"]], control="shuffle", objfun="de",
            meme_path="~/software/meme/bin/",
            outdir = dirname(stremeout_100bp_up))
}

print("Finished running streme for up- and down-regions")


#-------------------------------
# DREME
#-------------------------------

dremeout_100bp_down <- here("results/current/memes_bioc/dreme_100bp_down/dreme.xml")
if (!file.exists(dremeout_100bp_down)){
  runDreme(summit_flank_seq_bydirchange[["down"]], "shuffle",
           meme_path="~/software/meme/bin/",
           outdir = dirname(dremeout_100bp_down))
  }

dremeout_100bp_up <- here("results/current/memes_bioc/dreme_100bp_up/dreme.xml")
if (!file.exists(dremeout_100bp_up)){
  runDreme(summit_flank_seq_bydirchange[["up"]], "shuffle",
           meme_path="~/software/meme/bin/",
           outdir = dirname(dremeout_100bp_up))
  }

print("Done running DREME")

#-------------------------------
## run ame - discriminative mode
#-------------------------------

# enriched in upregulated with "down" as control
ame_discr_up <- here("results/current/memes_bioc/ame_discr_up/ame.tsv")
if (!file.exists( ame_discr_up )){
  runAme(summit_flank_seq_bydirchange, control = "down",
         meme_path=my_memepath,
         outdir=dirname(ame_discr_up))
}
# enriched in downregulated with "up" as control
ame_discr_down <- here("results/current/memes_bioc/ame_discr_down/ame.tsv")
if (!file.exists(ame_discr_down)){
  runAme(summit_flank_seq_bydirchange, control = "up",
         meme_path=my_memepath,
         outdir=dirname(ame_discr_down))
}

print("Finished running ame in discriminative mode for direction of expressionchange")


#-------------------------------
## all summits
#-------------------------------

## run streme to discover consensus motif
#-------------------------------
#print("Starting streme 1000bp")
#stremeout_1000bp <- here("results/current/memes_bioc/streme_1000bp/streme.xml")
#if (!file.exists( stremeout_1000bp )){
#  runStreme(summit_flank_1000bp_seq, control="shuffle",
#            meme_path="~/software/meme/bin/",
#            outdir = dirname(stremeout_1000bp) )
#}

print("Starting streme 100bp for all summits")
stremeout_100bp <- here("results/current/memes_bioc/streme_100bp/streme.xml")
if (!file.exists( stremeout_100bp )){
  runStreme(summit_flank_100bp_seq, control="shuffle",
            meme_path="~/software/meme/bin/",
            outdir = dirname(stremeout_100bp) )
}
print("Finished running streme for 100bp summit regions")

# Option objfun="cd" does not seem to get passed on to streme
#print("Starting streme 1000bp with central enrichment")
#stremeout_cd_1000bp <- here("results/current/memes_bioc/streme_cd_1000bp/streme.xml")
#if (!file.exists( stremeout_cd_1000bp )){
#  runStreme(summit_flank_1000bp_seq[1:100], objfun="cd", control=NA,
#            meme_path="~/software/meme/bin/",
#            outdir = dirname(stremeout_cd_1000bp) )
#}
#print("Finished running streme for 1000bp summit regions")

print("Analysis DONE")
