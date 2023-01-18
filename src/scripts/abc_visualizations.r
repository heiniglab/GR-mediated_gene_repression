suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))

option_list <- list(
  make_option(c( "--ABC_DexLPS_all"),
              type="character",
              help="Path to file with all ABC predictions passing 0.02 (in DexLPS condition)"),
  make_option(c("--ABC_LPS_all"),
              type="character",
              help="Path to file with all ABC predictions passing 0.02 (in LPS condition)"),
  make_option(c("--contrast_DexVSDexLPS"),
              type="character",
              help="Path to annotated tsv file of DeSeq2 contrast of DexLPS vs LPS"),
  make_option(c("--chipseq_summits"),
              type="character",
              help="Path to summit file of IDR peaks"),
  make_option(c("--igv"),
              type="character",
              help="Path to IGV snapshot of genomic locus (created manually)")
)

opt <- parse_args(OptionParser(option_list=option_list))

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages(library(dplyr, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(org.Mm.eg.db, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(DESeq2, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(ComplexHeatmap, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(TxDb.Mmusculus.UCSC.mm10.knownGene, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(topGO, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(ggExtra, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(plyr, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(here, warn.conflicts=F, quietly=T))

#set defaults for ggplot2 figures
theme_update(panel.background = element_rect(fill = "transparent", colour = NA),
             plot.background = element_rect(fill = "transparent", colour = NA),
             legend.background = element_rect(fill = "transparent", colour = NA),
             legend.key = element_rect(fill = "transparent", colour = NA),
             text=element_text(size=8, family = "ArialMT", colour="black"),
             title=element_text(size=10, family="ArialMT", colour="black"),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             axis.text = element_text(size=8, family="ArialMT", colour="black"),
             axis.line = element_line(colour="black"),
             axis.ticks = element_line(colour="black"))

#----------------------------------------------------------------------
##--NOTE: example code for additional plots can be found in abc_visualizations.Rmd
#----------------------------------------------------------------------

#----------------------------------------------------------------------
##----------- load data
#----------------------------------------------------------------------

chipseq_summits <- rtracklayer::import.bed(opt$chipseq_summits)

contrast_DexLPSvLPS <- read.delim(opt$contrast_DexVSDexLPS)
# rename first column for merge later
colnames(contrast_DexLPSvLPS)[1] <- "EnsemblID"

ABC_DexLPS_all <- read.delim(opt$ABC_DexLPS_all)
ABC_LPS_all <- read.delim(opt$ABC_LPS_all)
# remove version number from ensembl IDs
ABC_DexLPS_all$TargetGene = gsub("\\.[0-9]*$","",ABC_DexLPS_all$TargetGene)
ABC_LPS_all$TargetGene = gsub("\\.[0-9]*$","",ABC_LPS_all$TargetGene)


#----------------------------------------------------------------------
##----------- data exploration
#----------------------------------------------------------------------

# What do the distributions of ABC scores look like

ggplot()+
  geom_density(data=ABC_DexLPS_all %>% filter(!class=="promoter"), aes(x=ABC.Score, colour="DexLPS"))+
  geom_density(data=ABC_LPS_all %>% filter(!class=="promoter"), aes(x=ABC.Score, colour="LPS"))+
  scale_colour_manual("",
                      breaks = c("DexLPS", "LPS"),
                      values = c("#339966", "#0066CC"))

# How do the ABC scores of promoters compare to those of enhancers?
  
ggplot()+
geom_density(data=ABC_DexLPS_all %>% filter(class=="promoter"), aes(x=ABC.Score, colour="promoter"))+
geom_density(data=ABC_DexLPS_all %>% filter(class=="genic"), aes(x=ABC.Score, colour="genic"))+
geom_density(data=ABC_DexLPS_all %>% filter(class=="intergenic"), aes(x=ABC.Score, colour="intergenic"))+
scale_colour_manual("",
                    breaks = c("promoter", "genic", "intergenic"),
                    values = c("red", "green", "blue"))

# What distance is used to define sth as promoter?

ggplot()+
  geom_density(data=ABC_DexLPS_all, aes(x=distance, colour=class))+
  scale_x_continuous(trans = "log10")+
  scale_colour_manual("",
                      breaks = c("promoter", "genic", "intergenic"),
                      values = c("red", "green", "blue"))

ggplot()+
  geom_density(data=ABC_DexLPS_all %>% filter(class=="promoter"), aes(x=distance, colour=isSelfPromoter))+
  scale_x_continuous(trans = "log10")+
  scale_colour_manual("isSelfPromoter",
                      breaks = c("True", "False"),
                      values = c("red", "blue"))+
  labs(title="All promoter regions, split by whether they're promoters for their target gene")+
  theme(
    plot.title = element_text(size=12)
  )
#----------------------------------------------------------------------
# ------ ------  Number of enhancers per gene
#----------------------------------------------------------------------
# Are most genes regulated by a single enhancer region or by multiple? How many does the ABC model assign per gene?

plot_enhancers_per_gene <- function(){
  DexLPS_TargetGeneCounts <- plyr::count(ABC_DexLPS_all %>% filter(!class=="promoter"), vars="TargetGene")
  LPS_TargetGeneCounts <- plyr::count(ABC_LPS_all %>% filter(!class=="promoter"), vars="TargetGene")
  
  gg <- ggplot()+
    geom_histogram(data=DexLPS_TargetGeneCounts, aes(x=freq, fill="DexLPS", colour="DexLPS"), binwidth = 1,  alpha=0.2)+
    geom_histogram(data=LPS_TargetGeneCounts, aes(x=freq, fill="LPS", colour="LPS"), binwidth = 1, alpha=0.2)+
    scale_fill_manual("",
                      breaks = c("DexLPS", "LPS"),
                      values = c("#339966", "#0066CC")) +
    scale_colour_manual("",
                        breaks = c("DexLPS", "LPS"),
                        values = c("#339966", "#0066CC")) +
    guides(color = "none")+
    theme(legend.position = c(0.8,0.8))+
    labs(title=" ",
         x="# of enhancers per gene",
         y="counts")
  
  return(gg)
}

gg_enh_per_gene <- plot_enhancers_per_gene()
gg_enh_per_gene

print("We find an average of")
plyr::count(ABC_DexLPS_all %>% filter(!class=="promoter"), vars="TargetGene") %>% pull(freq) %>% mean()
print("enhancers per gene for the DexLPS condition and")
plyr::count(ABC_LPS_all %>% filter(!class=="promoter"), vars="TargetGene") %>% pull(freq) %>% mean()
print("in LPS.")

#----------------------------------------------------------------------
# ------ ------  Number of genes per enhancer
#----------------------------------------------------------------------
# ABC allows for multiple assignments (an individual enhancer assigned to more than one promoter)  

plot_genes_per_enhancer <- function(){
  
  DexLPS_EnhancerCounts <- plyr::count(ABC_DexLPS_all %>% filter(!class=="promoter"),vars="name")
  LPS_EnhancerCounts <- plyr::count(ABC_LPS_all %>% filter(!class=="promoter"),vars="name")
  
  gg_full <- ggplot()+
    geom_histogram(data=DexLPS_EnhancerCounts, aes(x=freq, fill="DexLPS",colour="DexLPS"), binwidth = 1, alpha=0.2)+
    geom_histogram(data=LPS_EnhancerCounts, aes(x=freq, fill="LPS",colour="LPS"), binwidth = 1, alpha=0.2)+
    scale_fill_manual("",
                      breaks = c("DexLPS", "LPS"),
                      values = c("#339966", "#0066CC")) +
    scale_colour_manual("",
                        breaks = c("DexLPS", "LPS"),
                        values = c("#339966", "#0066CC")) +
    guides(color = FALSE)+
    labs(title="Number of genes an individual enhancer is assigned to",
         x="# of genes per enhancer",
         y="counts")
  plot(gg_full)
  
  gg_zoom <- ggplot()+
    geom_histogram(data=DexLPS_EnhancerCounts, aes(x=freq, fill="DexLPS",colour="DexLPS"), binwidth = 1, alpha=0.2)+
    geom_histogram(data=LPS_EnhancerCounts, aes(x=freq, fill="LPS",colour="LPS"), binwidth = 1, alpha=0.2)+
    scale_fill_manual("",
                      breaks = c("DexLPS", "LPS"),
                      values = c("#339966", "#0066CC")) +
    scale_colour_manual("",
                        breaks = c("DexLPS", "LPS"),
                        values = c("#339966", "#0066CC")) +
    xlim(c(0,8))+
    guides(color = FALSE)+
    theme(legend.position = c(0.8,0.8))+
    labs(title=" ",
         x="# of genes per enhancer",
         y="counts")
  plot(gg_zoom)
  
  print("We find an average of")
  DexLPS_EnhancerCounts %>% pull(freq) %>% mean() %>% print()
  print("genes per enhancer in the DexLPS condition and")
  LPS_EnhancerCounts %>% pull(freq) %>% mean() %>% print()
  print("in LPS.")
  
  return(gg_zoom)
}

gg_genes_per_enhancer_zoom <- plot_genes_per_enhancer()


#----------------------------------------------------------------------
# ------ Comparing ABC Scores (including promoter regions) between the conditions
#----------------------------------------------------------------------

# merge regions from both conditions but also keep those that are only in one
# 
# add info no whether the target genes are DE or now and what direction they're changed

LPS_GR<- plyranges::as_granges(ABC_LPS_all, seqnames=chr)
DexLPS_GR <- plyranges::as_granges(ABC_DexLPS_all, seqnames=chr)

#----------------------------------------------------------------------
merge_LPSandDexLPS_ABCdata <- function(ABC_LPS_all_GR, ABC_DexLPS_all_GR){
  ABC_LPS_all_unique <- unique(ABC_LPS_all_GR)
  ABC_DexLPS_all_unique <- unique(ABC_DexLPS_all_GR)
  
  # and then figure out which ones overlap with one another
  overlap_ABC_all_unique <- ChIPpeakAnno::findOverlapsOfPeaks (ABC_LPS_all_unique,
                                                               ABC_DexLPS_all_unique,
                                                               connectedPeaks = "keepAll",
                                                               ignore.strand = TRUE)
  # we coerce the info of overlapping regions into a dataframe
  overlap_ABC_all_unique_df <- as.data.frame(overlap_ABC_all_unique$overlappingPeaks)
  
  # and clean up column names
  colnames(overlap_ABC_all_unique_df)<- gsub("ABC_LPS_all_unique...ABC_DexLPS_all_unique.","", colnames(overlap_ABC_all_unique_df))
  colnames(overlap_ABC_all_unique_df)<- gsub(".1", "", colnames(overlap_ABC_all_unique_df))
  colnames(overlap_ABC_all_unique_df)<- paste( c(rep("LPS_",27),rep("DexLPS_",27),rep("",2)), colnames(overlap_ABC_all_unique_df), sep="")
  
  # we use this overlap to define a new variable that holds the info of which regions of the conditions form a pair
  overlap_ABC_all_unique_df$pairID <- paste("pair", seq(1, nrow(overlap_ABC_all_unique_df)), sep="" )
  
  #-----------------for DexLPS-------------------
  # merge the pairID onto the original dataframes
  ABC_DexLPS_all_wpairID <- merge(ABC_DexLPS_all,
                                  overlap_ABC_all_unique_df[,c("DexLPS_name","pairID")],
                                  by.x="name",
                                  by.y="DexLPS_name",
                                  all.x=TRUE)
  # for those that don't correspond to a pair, we assign condition and name as pairID
  ABC_DexLPS_all_wpairID <- ABC_DexLPS_all_wpairID %>% 
    mutate(pairID = case_when(is.na(pairID) ~ paste("DexLPS",name, sep="_"),
                              TRUE ~ pairID) )
  # combination of pairID and gene they're assigned to forms the assignment variable.
  # we will use this later to plot values from the 2 conditions, that have the same assignment
  ABC_DexLPS_all_wpairID$assignment <- paste(ABC_DexLPS_all_wpairID$pairID,
                                             ABC_DexLPS_all_wpairID$TargetGene,
                                             sep="_")
  
  #-----------------for LPS-------------------
  # merge the pairID onto the original dataframes
  ABC_LPS_all_wpairID <- merge(ABC_LPS_all,
                               overlap_ABC_all_unique_df[,c("LPS_name","pairID")],
                               by.x="name",
                               by.y="LPS_name",
                               all.x=TRUE)
  # for those that don't correspond to a pair, we assign condition and name as pairID
  ABC_LPS_all_wpairID <- ABC_LPS_all_wpairID %>% 
    mutate(pairID = case_when(is.na(pairID) ~ paste("LPS",name, sep="_"),
                              TRUE ~ pairID) )
  
  ABC_LPS_all_wpairID$assignment <- paste(ABC_LPS_all_wpairID$pairID,
                                          ABC_LPS_all_wpairID$TargetGene,
                                          sep="_")
  
  merged_conditions <- merge(ABC_DexLPS_all_wpairID,
                                 ABC_LPS_all_wpairID,
                                 by="assignment", 
                                 all=TRUE,
                                 suffixes = c(".DexLPS", ".LPS"))
  merged_conditions <- merged_conditions %>%
    mutate(TargetGene = case_when(TargetGene.DexLPS==TargetGene.LPS ~ TargetGene.DexLPS,
                                  is.na(TargetGene.DexLPS) ~ TargetGene.LPS,
                                  is.na(TargetGene.LPS) ~ TargetGene.DexLPS,
                                  TRUE ~ NA_character_)
    )
  
  return(merged_conditions)
}

merged_conditions_ABC <- merge_LPSandDexLPS_ABCdata(LPS_GR,DexLPS_GR)

#----------------------------------------------------------------------
# ------ merge on gene expression results
#----------------------------------------------------------------------

# merge the RNAseq results to use as colour later
merged_conditions_ABC <- merge(merged_conditions_ABC, contrast_DexLPSvLPS, 
                               by.x="TargetGene", by.y="EnsemblID")
merged_conditions_ABC <- merged_conditions_ABC %>%
  mutate(change=case_when(padj<0.05 & log2FoldChange>0.58 ~ "up",
                          padj<0.05 & log2FoldChange<(-0.58) ~ "down",
                          TRUE ~ "no change"))

# replace the na with 0
merged_conditions_ABCnomissing <- merged_conditions_ABC %>% 
  mutate(ABC.Score.DexLPS = case_when(is.na(ABC.Score.DexLPS) ~ 0,
                                      TRUE ~ ABC.Score.DexLPS) ) %>%
  mutate(ABC.Score.LPS = case_when(is.na(ABC.Score.LPS) ~ 0,
                                   TRUE ~ ABC.Score.LPS) )

cor_ABCscores <- cor(merged_conditions_ABCnomissing$ABC.Score.DexLPS,
                     merged_conditions_ABCnomissing$ABC.Score.LPS)

gg_merged_conditions_ABCnomissing <- ggplot(merged_conditions_ABCnomissing)+
  geom_point(aes(x=ABC.Score.DexLPS, y=ABC.Score.LPS, fill=log2FoldChange), alpha=0.2, size=1, stroke = 0, shape=21)+
  scale_fill_gradient(low="blue", high="red")+
  ylim(c(0,0.8))+
  xlim(c(0,0.8))+
  geom_text(x=0.2, y=0.75, label=paste("r = ",format(cor_ABCscores, digits = 2)))+
  labs(fill="log2FC", x="ABC score Dex+LPS", y="ABC score LPS")+
  theme(legend.position = "top")

gg_merged_conditions_ABCnomissing

#----------------------------------------------------------------------
# ------ plot foldchange in ABC score vs foldchange in expression
#----------------------------------------------------------------------

cor_ABCdiff_log2FC <- with(merged_conditions_ABCnomissing %>% filter(change!="no change"),
                           cor((ABC.Score.DexLPS - ABC.Score.LPS),log2FoldChange))

gg_ABCdiff_log2FC <- ggplot(merged_conditions_ABCnomissing %>% filter(change!="no change"), aes(x=(ABC.Score.DexLPS - ABC.Score.LPS), y=log2FoldChange)) +
  geom_hex(bins=70)+
  viridis::scale_fill_viridis()+
  geom_smooth(method = "lm", colour="black")+
  geom_text(x=0.6, y=5, label=paste("r = ",format(cor_ABCdiff_log2FC, digits = 2)))+
  labs(x="ABC score Dex+LPS - ABC score LPS", y="expression log2FC(Dex+LPS / LPS)")+
  theme(legend.position = "top")



#----------------------------------------------------------------------
# ------ ABC score by peak presence
#----------------------------------------------------------------------

plot_ABC_with_marginals <- function(whatchange){
  gg <- ggplot(data= merged_conditions_ABCnomissing %>% filter(change==whatchange) , 
               aes(x=ABC.Score.DexLPS, y=ABC.Score.LPS ) )+
    geom_point( aes( fill=whatchange) , size=1, alpha=1,stroke = 0, shape=21) +
    geom_abline(slope = 1)+
    ylim(c(0,0.8))+
    xlim(c(0,0.8))+
    labs(x="ABC score Dex+LPS", y="ABC score LPS")+
    scale_fill_manual(
      breaks=c("up","no change","down"),
      values=c("red","black","blue") )+
    theme(legend.position = "none")
  
  gg_marg <- ggExtra::ggMarginal(gg, type="histogram", size=2,
                                 xparams = list(bins=85), yparams = list(bins=85))
  
  return(gg_marg)
}

ggplot_ABC_with_marginals_up <- plot_ABC_with_marginals("up")
ggplot_ABC_with_marginals_up
ggplot_ABC_with_marginals_down <-plot_ABC_with_marginals("down")
ggplot_ABC_with_marginals_down

#----------------------------------------------------------------------
# ------ show marginals as bar plot
#----------------------------------------------------------------------
gg_marginals_dexlps <- 
  ggplot()+
  geom_histogram(data=merged_conditions_ABCnomissing %>% filter(change=="up"), aes(x=ABC.Score.DexLPS, fill="up"), alpha=0.2, bins=85)+
  geom_histogram(data=merged_conditions_ABCnomissing %>% filter(change=="down"), aes(x=ABC.Score.DexLPS, fill="down"), alpha=0.2, bins=85)+
  scale_fill_manual(
    breaks=c("up","down"),
    values=c("red","blue") )+
  xlim(c(NA,0.3))+
  labs(fill="gene expression change",x="ABC score DexLPS", y="ABC score LPS")+
  theme(legend.position = c(0.7,0.7))
gg_marginals_lps <- 
  ggplot()+
  geom_histogram(data=merged_conditions_ABCnomissing %>% filter(change=="up"), aes(x=ABC.Score.LPS, fill="up"), alpha=0.2, bins=85)+
  geom_histogram(data=merged_conditions_ABCnomissing %>% filter(change=="down"), aes(x=ABC.Score.LPS, fill="down"), alpha=0.2, bins=85)+
  scale_fill_manual(
    breaks=c("up","down"),
    values=c("red","blue") )+ 
  xlim(c(NA,0.3))+
  labs(fill="gene expression change",x="ABC score DexLPS", y="ABC score LPS")+
  theme(legend.position = c(0.7,0.7))

#----------------------------------------------------------------------
# ------ ABC score by peak presence
#----------------------------------------------------------------------
# Is there a difference in ABC score between the enhancers overlapping with a GR summit vs the ones that don't?


DexLPS_GR_summits_overlap <- IRanges::findOverlaps(DexLPS_GR,chipseq_summits) #(query,subject)
DexLPS_GR$haspeakID <- FALSE
DexLPS_GR$haspeakID[queryHits(DexLPS_GR_summits_overlap)] <- TRUE


plot_ABCscore_ofDEgenes_byGRoverlap <- function(){
  DE_ENSEMBL <- contrast_DexLPSvLPS %>%  filter(padj<0.05 & abs(log2FoldChange)>0.58 ) %>% pull(EnsemblID)
  
  DexLPS_enhancers_DEsubset <- 
    as.data.frame(DexLPS_GR) %>% 
    filter(TargetGene %in% DE_ENSEMBL) %>%
    filter(!class=="promoter")
  
  DexLPS_enhancers_DEsubset %>% filter(haspeakID==TRUE) %>% dplyr::pull(ABC.Score) %>% mean() %>% print()
  DexLPS_enhancers_DEsubset %>% filter(haspeakID==FALSE) %>% dplyr::pull(ABC.Score) %>% mean() %>% print()
  
  # show distribution for the two groups
  
  gg1 <- ggplot()+
    geom_density(data=DexLPS_enhancers_DEsubset %>% filter(haspeakID==TRUE), 
                 aes(x=ABC.Score, colour="haspeakID", fill="haspeakID"), alpha=0.3)+
    geom_density(data=DexLPS_enhancers_DEsubset %>% filter(haspeakID==FALSE), 
                 aes(x=ABC.Score, colour="nopeakID", fill="nopeakID"), alpha=0.3)+
    scale_x_continuous(trans = "log2")+
    scale_colour_manual("",
                        labels = c("has GR peak", "no GR peak"),
                        breaks = c("haspeakID", "nopeakID"),
                        values = c("darkmagenta", "cadetblue"))+
    scale_fill_manual("",
                      labels = c("has GR peak", "no GR peak"),
                      breaks = c("haspeakID", "nopeakID"),
                      values = c("darkmagenta", "cadetblue"))+
    guides(colour="none")+
    #xlim(c(-3,0))+
    labs(title=" ",
         x="ABC score Dex+LPS",
         y="density")+
    theme(
      axis.text = element_text(size=12, family="ArialMT", colour="black"),
      axis.title = element_text(size=12, family="ArialMT", colour="black"),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
      legend.position = c(0.8,0.8)
    )
  plot(gg1)
  
  gg2 <- DexLPS_enhancers_DEsubset %>%
    ggplot( aes(x=hic_contact_pl_scaled_adj, y=log2(activity_base))) +
    geom_hex()+
    facet_wrap(~haspeakID, labeller="label_both")+
    viridis::scale_fill_viridis()+
    labs(title = " ", x="Hi-C contact", y="log2(base activity)")
  
  plot(gg2)
  all_plots=list()
  all_plots[[1]] <- gg1
  all_plots[[2]] <- gg2
  return(all_plots)
}

gg_ABCscore_ofDEgenes_byGRoverlap <- plot_ABCscore_ofDEgenes_byGRoverlap()


#----------------------------------------------------------------------
# ------ load IGV plot
#----------------------------------------------------------------------
igv <- png::readPNG( opt$igv )

gg_igv <- ggplot() + 
  ggpubr::background_image(igv) +
  # so it doesn't get squished
  #coord_fixed()+
  # This ensures that the image leaves some space at the edges
  theme(plot.margin = margin(t=0, l=0, r=0, b=0, unit = "cm"),
        axis.line = element_blank())


#----------------------------------------------------------------------
# ------ put figure panel together
#----------------------------------------------------------------------

# save giant files separately
ggsave(here("./results/current/Figures/Figure_abcresults_3B.bmp"), gg_merged_conditions_ABCnomissing,
       width=75, height=100, units="mm",
       bg="white")
ggsave(here("./results/current/Figures/Figure_abcresults_3B.png"), gg_merged_conditions_ABCnomissing,
       width=75, height=100, units="mm",
       bg="white")

gg_placeholder <- ggplot() + 
  theme(plot.margin = margin(t=0, l=0, r=0, b=0, unit = "cm"),
        axis.line = element_blank())

# first row
gg_r1c1 <-  ggpubr::ggarrange(gg_enh_per_gene, gg_genes_per_enhancer_zoom,
                            labels = c("A", NA),
                            ncol = 1, nrow = 2, heights = c(1,1))

gg_r1c2 <-  ggpubr::ggarrange(gg_placeholder, ggplot_ABC_with_marginals_up, ggplot_ABC_with_marginals_down, 
                            labels = c("B", "C", NA),
                            ncol = 3, nrow=1, widths=c(1,1,1))

gg_r1 <- ggpubr::ggarrange(gg_r1c1, gg_r1c2, 
                           labels = c(NA, NA),
                           ncol = 2, nrow=1, widths=c(1,3))

# second row
gg_r2 <-  ggpubr::ggarrange(gg_ABCdiff_log2FC, gg_ABCscore_ofDEgenes_byGRoverlap[[1]], gg_ABCscore_ofDEgenes_byGRoverlap[[2]], 
                            labels = c("D","E", "F"),
                            ncol = 3, nrow=1, widths=c(1,1,1.5))


gg_tophalf <- ggpubr::ggarrange(gg_r1, gg_r2, 
                                nrow=2, heights =c(1,1))
gg_tophalf

full_panel <- ggpubr::ggarrange(gg_tophalf, gg_igv, 
                                labels=c(NA,"G"),
                                nrow=2, heights = c(1.8,1))

#full_panel
ggsave(here("./results/current/Figures/Figure_abcresults.png"), full_panel,
       width=300, height=300, units="mm",
       bg="white")
ggsave(here("./results/current/Figures/Figure_abcresults.pdf"), full_panel,
       width=300, height=300, units="mm",
       bg="white")

