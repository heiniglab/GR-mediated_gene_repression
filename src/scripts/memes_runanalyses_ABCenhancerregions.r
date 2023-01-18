suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))

option_list <- list(
  make_option(c("--ABC_all"), 
              type="character",
              help="path to abc results of dexlps condition"),
  make_option(c("--memedb_expressed"),
              type="character",
              help="Path to memedb file filtered for motifs where TFs are expressed in 4sU"),
  make_option(c("--output"),
              type="character",
              help="fimo results file")
  )

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
ABC_all <- read.delim(opt$ABC_all) %>% plyranges::as_granges(., seqnames=chr)

# no need to run fimo a bunch of times on the same enhancers regions, just because they are listed more than once (with different ABCscores)
ABC_unique <- unique(ABC_all)

meme_db_expressed <- readRDS(opt$memedb_expressed)

# to_list() converts the database back from data.frame format to a standard `universalmotif` object.
options(meme_db = to_list(meme_db_expressed, extrainfo = FALSE))

# where is meme installed 
my_memepath="~/software/meme/bin/"
check_meme_install(meme_path=my_memepath)

#-------------------------------
## get sequences
#-------------------------------

enhancer_seq <- ABC_unique %>%
  get_sequence(mm.genome)

#-------------------------------
## run fimo
#-------------------------------

# conda activate py_3
# perlbrew use perl-5.34.0
# nohup Rscript memes_runanalyses_ABCenhancerregions.r & (from within the script directory)

fimo_results <-
  runFimo(enhancer_seq,
          meme_db_expressed,
          meme_path=my_memepath)

print("Finished running fimo")

saveRDS (fimo_results, opt$output)

print("Done saving results")
