suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))

option_list <- list(
  make_option(c("-s", "--peaksummits_annot_resDexLPSvLPS"),
              type="character",
              help="Path to annotated peaksummits"),
  make_option(c("-m", "--fimo_featurematrix"),
              type="character",
              help="Path to matrix with motifcounts for peakregions"),
  make_option(c("-o", "--out_assign"),
              type="character",
              help="Path of outfile with peakregion to gene assignments")
)

suppressPackageStartupMessages(library(dplyr, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts=F, quietly=T))

opt <- parse_args(OptionParser(option_list=option_list))


#set defaults for ggplot2 figures
theme_update(panel.background = element_rect(fill = "transparent", colour = NA),
             plot.background = element_rect(fill = "transparent", colour = NA),
             legend.background = element_rect(fill = "transparent", colour = NA),
             legend.key = element_rect(fill = "transparent", colour = NA),
             text=element_text(size=6, family = "ArialMT", colour="black"),
             title=element_text(size=8, family="ArialMT", colour="black"),
             panel.grid.major = element_line(colour="grey", size=0.2),
             panel.grid.minor = element_blank(),
             axis.text = element_text(size=6, family="ArialMT", colour="black"),
             axis.line = element_line(colour="black"),
             axis.ticks = element_line(colour="black"),
             legend.key.size = unit(6, 'points'), #change legend key size
             legend.key.height = unit(6, 'points'), #change legend key height
             legend.key.width = unit(6, 'points'), #change legend key width
             legend.text = element_text(size=6, family="ArialMT", colour="black")
)

# Load our original annotation (done by proximity)
summits_annot2ALLgenes_resDexLPSvLPS <- read.delim(opt$peaksummits_annot_resDexLPSvLPS)

# Load motifdata, to see what peak summits are within the same peakregion
motifdata <- read.delim(opt$fimo_featurematrix)

### Encode peak - gene assignment as separate variable

# For the ABC annotation, a single enhancer can map to multiple genes, for the proximity based one only to one.  
# In order to compare the results, we make the assignment a new variable

# make a new variable out of peakID - gene combo?
summits_annot2ALLgenes_resDexLPSvLPS <- summits_annot2ALLgenes_resDexLPSvLPS %>%
  mutate(assignment = paste(name,ENSEMBL, sep="_"))

# dataframe for the proximity based annotations
assignment_prox <- data.frame(
  peakID = character(),
  anno = character()
)

# motifdata is the data holding features for each of the input peakregions
for (peakregion_name in motifdata$name) {
  # some peakregions overlap multiple peaksummits (peaks that got merged cuz they overlap) we can split them based on "|"
  peakIDs <- strsplit(peakregion_name,"[|]")[[1]] %>% as.vector()
  
  #get the corresponding proximity-based annotaion
  #--------------------------------------------------
  prox_anno <- summits_annot2ALLgenes_resDexLPSvLPS %>% 
    filter(name %in% peakIDs) %>% 
    pull(ENSEMBL)
  
  # check whether the annotations are the same by confirming that the average # of annotations matches the unique ones
  avg_prox_anno_no <- length(prox_anno) / length(peakIDs)
  if ( avg_prox_anno_no != length(unique(prox_anno)) ){
    print(paste(paste(peakIDs, collapse=" and "),"annotate to different genes" ))
  }
  
  assignment_prox_new <- data.frame(peakID=rep(peakregion_name,length(unique(prox_anno))),
                                    anno = unique(prox_anno))
  assignment_prox <- rbind(assignment_prox,assignment_prox_new)
  
}

# There's about a dozen cases in which peaks are merged, but their summits didn't annotate to the same gene. 
# ATM we just combine the annotations for both if their peaks got merged. 

write.table(assignment_prox, file=opt$out_assign, sep="\t", quote = FALSE, row.names = FALSE)
