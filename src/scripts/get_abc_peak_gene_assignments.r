suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))

option_list <- list(
  make_option(c("-a", "--ABCscores_ftr0.02"),
            type="character",
            help="Path to filtered ABC scores"),
  make_option(c("-s", "--peaksummits"),
            type="character",
            help="Path topeaksummits of reliable peaks"),
  make_option(c("-m", "--fimo_featurematrix"),
            type="character",
            help="Path to matrix with motifcounts for peakregions"),
  make_option(c("-o", "--out_assignment_peakID"),
            type="character",
            help="Path of outfile with peakID to gene assignments"),
  make_option(c("-r", "--out_assignment_peakregion"),
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

# Read in ABC scores that passed threshold of 0.02
ABCscores_ftr0.02 <- read.delim(opt$ABCscores_ftr0.02)

# Load the "reliable peak" coordinates.  
summits <- rtracklayer::import.bed(opt$peaksummits)

# Load motifdata, to see what peak summits are within the same peakregion
motifdata <- read.delim(opt$fimo_featurematrix)

#-----------------------------------------------------------------------------------
# Get ABC-based peak-gene assignments per condition
#-----------------------------------------------------------------------------------

# In order to know what peakID gets annotated to what gene, we need to annotate enhancers to peakIDs first by seeing what enhancer the peak summits fall into.  
# Once we have this assignment, we can check whether the ABC enhancer region and the peak summit are annotated to the same gene.

### 1. See what enhancer the summit falls into.

# we can't keep duplicated ranges anyways, so we can toss the extra columns and then simply use the name to match the coordinates back to their genes
uniqueranges_ABCresults <- unique(ABCscores_ftr0.02[,c(1:4)])

uniqueranges_ABCresults <- GenomicRanges::makeGRangesFromDataFrame(uniqueranges_ABCresults,
                                                             keep.extra.columns = TRUE)

summits_uniqueranges_overlap <- ChIPpeakAnno::findOverlapsOfPeaks (summits,
                                                                uniqueranges_ABCresults,
                                                                connectedPeaks = "keepAll",
                                                                ignore.strand = TRUE)

### 2. Find ABC model derived assignments for our peaks  

# Grabbing the assignment of summit name to enhancer name and through that assign the summit to the respective ABC derived target gene.  
# For an enhancer mapping to multiple genes, the summits should get multiple rows: each for one gene assignment

# transform into dataframe and clean up the column names
summits_uniqueranges_overlap_df <- summits_uniqueranges_overlap$overlappingPeaks %>% 
  as.data.frame()
colnames(summits_uniqueranges_overlap_df)<- gsub("summits...uniqueranges_ABCresults.","",colnames(summits_uniqueranges_overlap_df))

# peaks1 corresponds to the summitIDs we created earlier, the name belongs to the enhancers
summit_peak_assignment <- summits_uniqueranges_overlap_df %>% dplyr::select("name","start","end","name.1")
colnames(summit_peak_assignment) <- c("peakID","summit_start","summit_end","enhancerID")

# now we can merge the summitIDs with the original ABC results based on the name of the enhancers
ABCscores_ftr0.02_assignedsummits <- merge(summit_peak_assignment, ABCscores_ftr0.02,
                                                        by.x="enhancerID",by.y="name",
                                                        all=TRUE)

### 3. Encode peak - gene assignment as separate variable

# For the ABC annotation, a single enhancer can map to multiple genes, for the proximity based one only to one.  
# In order to compare the results, we make the assignment a new variable

ABCscores_ftr0.02_assignedsummits_wpeakID <- ABCscores_ftr0.02_assignedsummits %>%
  filter(!is.na(peakID)) %>%
  mutate(assignment = paste(peakID,TargetGene, sep="_"))

### 4. Get annotation on per regions level (instead of per peaksummit)
# In order to merge our annotations with the featurematrix later, we need to accomodate those peakregions that contain more than one peaksummit
# I iterate through the peakregions and get annotations on a per region level (instead of on a peaklevel)
# If peak1 and peak2 are so close, that they got merged, they will map to largely the same enhancer region and thereby get the same gene annotations  
# In case they merge to different genes, we simply merge their annotations
# When we aggregate the feature later, we still wanna be able to use the ABC scores. Since there's multiple peak summits mapping to the same ABC enhancer regions we keep info of the ABC scores and can use it to weight the features later.

# dataframe for the ABC derived annotations including promoter regions
assignment_abc_wprom <- data.frame(
  peakID = character(),
  anno = character(),
  abcscore = character(),
  abcnumerator= character()
)

# motifdata is the data holding features for each of the input peakregions
for (peakregion_name in motifdata$name) {
  # some peakregions overlap multiple peaksummits (peaks that got merged cuz they overlap) we can split them based on "|"
  peakIDs <- strsplit(peakregion_name,"[|]")[[1]] %>% as.vector()

  #  get the ABC derived annotations (including promoters)
  #---------------------------------------
  abc_prom_anno <- ABCscores_ftr0.02_assignedsummits_wpeakID %>%
    filter(peakID %in% peakIDs) %>%
    dplyr::select("peakID","TargetGene","ABC.Score","ABC.Score.Numerator","class","isSelfPromoter") %>%
    distinct()
  
  assignment_abc_prom_new <- data.frame(peakID=rep(peakregion_name,nrow(abc_prom_anno)),
                                        anno = abc_prom_anno$TargetGene,
                                        abcscore = abc_prom_anno$ABC.Score,
                                        abcnumerator = abc_prom_anno$ABC.Score.Numerator,
                                        class  = abc_prom_anno$class,
                                        isSelfPromoter = abc_prom_anno$isSelfPromoter )
  assignment_abc_wprom <- rbind(assignment_abc_wprom,assignment_abc_prom_new)
  
  # check whether the annotations are the same by confirming that the average # of annotations matches the unique ones
  avg_abc_prom_anno_no <- nrow(abc_prom_anno) / length(peakIDs)
  if ( avg_abc_prom_anno_no != length(unique(abc_prom_anno$TargetGene)) ){
    print(paste(paste(peakIDs, collapse=" and "),"annotate to different genes" ))
  }
  
}


# There's about a dozen cases in which peaks are merged, but their summits didn't annotate to the same gene. 
# ATM we just combine the annotations for both if their peaks got merged. 
write.table(ABCscores_ftr0.02_assignedsummits,
            file=opt$out_assignment_peakID, 
            sep="\t", quote = FALSE, row.names = FALSE)

write.table(assignment_abc_wprom, 
            file=opt$out_assignment_peakregion, 
            sep="\t", quote = FALSE, row.names = FALSE)

