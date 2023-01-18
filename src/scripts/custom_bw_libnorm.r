suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))

option_list <- list(
  make_option(c( "--bw_DexLPS_h3k27ac"),
              type="character",
              help="Path to H3K27ac bigwig file from DexLPS condition"),
  make_option(c("--bw_LPS_h3k27ac"),
              type="character",
              help="Path to H3K27ac bigwig file from LPS condition"),
  make_option(c("--counts_h3k27ac"),
              type="character",
              help="Path to file with adjusted libsize (number of reads in bam file overlappign joint peak universe)"),
  make_option(c( "--bw_DexLPS_atac"),
              type="character",
              help="Path to ATAC bigwig file from DexLPS condition"),
  make_option(c("--bw_LPS_atac"),
              type="character",
              help="Path to ATAC bigwig file from LPS condition"),
  make_option(c("--counts_atac"),
              type="character",
              help="Path to file with adjusted libsize (number of reads in bam file overlappign joint peak universe)"),
  make_option(c("--gtf"),
              type="character",
              help="Path to gtf file of genomic reference")
  )

opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(library(here, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(dplyr, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(plyranges, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(rtracklayer, warn.conflicts=F, quietly=T))

#-----------------------------------------
#------- H3K27ac
#-----------------------------------------

bw_DexLPS_h3k27ac <- import.bw( opt$bw_DexLPS_h3k27ac )
bw_LPS_h3k27ac <- import.bw( opt$bw_LPS_h3k27ac )
counts_h3k27ac <- read.table(opt$counts_h3k27ac)

# normalization is important in this case as seen by the ratio
as.numeric(counts_h3k27ac[2,1]) / as.numeric(counts_h3k27ac[4,1])

# scale by "adjusted" library size (number of reads overlapping peak universe)
score(bw_DexLPS_h3k27ac) <- score(bw_DexLPS_h3k27ac) / ( as.numeric(counts_h3k27ac[2,1]) / 10^6)
score(bw_LPS_h3k27ac) <- score(bw_LPS_h3k27ac) / ( as.numeric(counts_h3k27ac[4,1]) / 10^6)

# reexport the normalized tracks
export.bw(bw_DexLPS_h3k27ac,
          here("./results/current/ChIP/H3K27ac/bw/DexLPS_histone_H3K27ac_PE_merged_GRCm38_libnorm.bw"))

export.bw(bw_LPS_h3k27ac,
          here("./results/current/ChIP/H3K27ac/bw/LPS_histone_H3K27ac_PE_merged_GRCm38_libnorm.bw"))

rm(list=c("bw_DexLPS_h3k27ac", "bw_LPS_h3k27ac"))

#-----------------------------------------
#------- ATAC
#-----------------------------------------

bw_DexLPS_atac <- import.bw( opt$bw_DexLPS_atac )
bw_LPS_atac <- import.bw( opt$bw_LPS_atac )

counts_atac <- read.table(opt$counts_atac)

# normalization is NOT so important in this case as seen by the ratio
as.numeric(counts_atac[2,1]) / as.numeric(counts_atac[4,1])

# scale by "adjusted" library size (number of reads overlapping peak universe)
score(bw_DexLPS_atac) <- score(bw_DexLPS_atac) / ( as.numeric(counts_atac[2,1]) / 10^6)
score(bw_LPS_atac) <- score(bw_LPS_atac) / ( as.numeric(counts_atac[4,1]) / 10^6)

# reexport the normalized tracks
export.bw(bw_DexLPS_atac,
          here("./results/current/atacseq/bw/merged_DexLPS_GRCm38_libnorm.bw"))

export.bw(bw_LPS_atac,
          here("./results/current/atacseq/bw/merged_LPS_GRCm38_libnorm.bw"))

#-----------------------------------------
#------- reference gtf file for IGV
#-----------------------------------------
#import gtf file with rtracklayer
gtf <- rtracklayer::import(opt$gtf)
gtf_nogene <- gtf %>% filter(type!="gene")
rtracklayer::export(gtf_nogene,
                    here("./data/current/Mus_musculus.GRCm38.100_ftrnogene.gtf"))
rtracklayer::export.gff3(gtf_nogene,
                         here("./data/current/Mus_musculus.GRCm38.100_ftrnogene.gff3"))
