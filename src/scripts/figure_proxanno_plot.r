suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))

option_list <- list(
  make_option(c( "--summitAnno_expr"),
              type="character",
              help="Path to file with GR idr summits annotated to expressed genes (as anno object)"),
  make_option(c("--summitAnno_df_expr"),
              type="character",
              help="Path to file with GR idr summits annotated to expressed genes (as dataframe)"),
  make_option(c("--permtest_res"),
              type="character",
              help="Path to rds file with results of permutation test of group differences"),
  make_option(c("--fimo_results"),
              type="character",
              help="Path to rds file with fimo hits within summitregions (1000bp)"),
  make_option(c("--chipseq_summit_granges"),
              type="character",
              help="Path to summit file of IDR peaks"),
  make_option(c("--deeptools"),
              type="character",
              help="Path to png of deeptoolsheatmap split by activating vs repressing peaks"),
  make_option(c("--streme_100bp_up"),
              type="character",
              help="Path to streme xml file for 100bp regions around activating peak regions"),
  make_option(c( "--streme_100bp_down"),
              type="character",
              help="Path to streme xml file for 100bp regions around repressing peak regions"))

opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(library(ggplot2, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(here, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(memes, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(universalmotif, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(dplyr, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(png, warn.conflicts=F, quietly=T))

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
             axis.ticks = element_line(colour="black"),
             legend.key.size = unit(8, 'points'), #change legend key size
             legend.key.height = unit(8, 'points'), #change legend key height
             legend.key.width = unit(8, 'points'), #change legend key width
             legend.text = element_text(size=8, family="ArialMT", colour="black")
)

#-------------------------------
## 0. load data
#-------------------------------

summitAnno_expr <- readRDS( opt$summitAnno_expr)
summitAnno_df_expr <- readRDS( opt$summitAnno_df_expr )
permtest_res <- readRDS( opt$permtest_res )
fimo_results <- readRDS( opt$fimo_results )
ChIPseq_summit_Granges <- readRDS( opt$chipseq_summit_granges )

# Take 100bp windows around ChIP-seq summits
summit_flank_100bp <- ChIPseq_summit_Granges %>% 
  plyranges::anchor_center() %>% 
  plyranges::mutate(width = 100) 

# Take 100bp windows around ChIP-seq summits
summit_flank_1000bp <- ChIPseq_summit_Granges %>% 
  plyranges::anchor_center() %>% 
  plyranges::mutate(width = 1000) 


#---------------------------------
#### Peak distance from TSS
#---------------------------------

#Turning things around, let's check 
#* how far the closest peak is from the TSS of up-vs downregulated gene  
#* how far we need to go from the TSS to have 2 or 3 peaks mapping to the gene

distbygene <-  summitAnno_df_expr  %>% 
  group_by(SYMBOL, change) %>%
  summarise(min_dist=min(abs(distanceToTSS)), 
            mean_dist=mean(abs(distanceToTSS))) %>%
  ungroup() %>%
  mutate(logmindist=log2(min_dist+1))

ggplot( ) +
  geom_histogram(data=distbygene %>% filter(change=="up"), aes(x=min_dist, fill="up"), alpha=.3, binwidth = 5000) +
  geom_histogram(data=distbygene %>% filter(change=="down"), aes(x=min_dist, fill="down"), alpha=.3,  binwidth = 5000) +
  expand_limits(x=c(0,400000))+
  labs(title="Dist. to nearest peak: split by direction of expression change")+
  scale_x_continuous(breaks = seq(0, 400000, by = 100000),
                     labels = paste(seq(0, 400000, by = 100000) / 1000,"kb"))+
  scale_fill_manual(name = "", values = c( "up" = "red", "down" = "blue"))

# Zoom into reasonable region

gg_distbychange <- ggplot( ) +
  geom_histogram(data=distbygene %>% filter(change=="ns"), aes(x=logmindist,y=..density.., fill="ns"), binwidth=0.5, alpha=.2) +
  geom_histogram(data=distbygene %>% filter(change=="up"), aes(x=logmindist,y=..density.., fill="up"), binwidth=0.5, alpha=.2) +
  geom_histogram(data=distbygene %>% filter(change=="down"), aes(x=logmindist,y=..density.., fill="down"), binwidth=0.5, alpha=.2) +
  
  geom_density(data=distbygene %>% filter(change=="up"), aes(x=logmindist, colour="up"), alpha=.3, show.legend = FALSE) +
  geom_density(data=distbygene %>% filter(change=="down"), aes(x=logmindist, colour="down"), alpha=.3, show.legend = FALSE) +
  geom_density(data=distbygene %>% filter(change=="ns"), aes(x=logmindist, colour="ns"), alpha=.3, show.legend = FALSE) +
  
  geom_vline(xintercept = log2(30000), size=1, colour="black", linetype=2)+
  scale_x_continuous(breaks = c( log2(0+1), log2(1000+1), log2(5000+1), log2(10000+1), log2(30000+1),log2(100000+1) ),
                     labels = paste(c(0,1000,5000,10000,30000,100000) / 1000,"kb"))+
  scale_fill_manual(name = "", values = c( "up" = "red", "down" = "blue", "ns" = "black"))+
  scale_colour_manual(name = "", values = c( "up" = "red", "down" = "blue", "ns" = "black"))+
  guides(colour="none")+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        legend.position = c(0.15,0.8))+
  labs(x="Distance to nearest peak - genecentric")

gg_distbychange

#---------------------------------
# --- load deeptools figure
#---------------------------------

deeptools <- png::readPNG( opt$deeptools )

gg_deeptools <- ggplot() + 
  ggpubr::background_image(deeptools) +
  # so it doesn't get squished
  #coord_fixed()+
  # This ensures that the image leaves some space at the edges
  theme(plot.margin = margin(t=0, l=0, r=0, b=0, unit = "cm"),
        axis.line = element_blank())


#--------------------------------------
#- permutations
#--------------------------------------
groupdiff <- as.numeric(permtest_res[[3]])

gg_perm_dist <- ggplot()+
  geom_histogram(aes(x=permtest_res[[2]]), bins=100)+
  labs(title="",
       x="Permutation based \n group differences")+
  xlim(-10000,10000)+
  geom_vline( xintercept = groupdiff, size=1, colour="red", linetype=1)

#--------------------------------------
#- memes
#--------------------------------------

streme_up_results <- memes::importStremeXML( opt$streme_100bp_up )
streme_up_results <- streme_up_results %>% 
  mutate(name = paste(consensus, pval))
                                     
gg_streme_up <-
  streme_up_results[1:5,] %>% 
  to_list() %>% 
  view_motifs(text.size = 7,
              tryRC = TRUE)
  #theme(plot.margin = margin(0,0,2.4,0, "cm"))
gg_streme_up

streme_down_results <- memes::importStremeXML( opt$streme_100bp_down )
streme_down_results <- streme_down_results %>% 
  mutate(name = paste(consensus, pval))

gg_streme_down <-
  streme_down_results[1:5,] %>% 
  to_list() %>% 
  view_motifs(text.size = 7,
              tryRC = TRUE)
gg_streme_down

#---------------------------------
# --- fimo chi-square
#---------------------------------

make_fimo_chisquare_plot <- function(summit_flank, fimo_results, seqwidth){
  input_intersect_hits <- summit_flank %>% 
    plyranges::join_overlap_intersect(fimo_results)
  
  fimo_counts <- 
    as.data.frame(input_intersect_hits) %>% 
    group_by(motif_alt_id,directionchange) %>% 
    summarize(count=n())
  
  fimo_counts <- fimo_counts %>% tidyr::pivot_wider(
    names_from = directionchange, 
    values_from = count)
  
  tbl_directionchange <- table(summit_flank$directionchange)
  setsize <- tbl_directionchange[["down"]] + tbl_directionchange[["up"]]
  
  fimo_counts <- fimo_counts %>%
    na.omit()
  
  fimo_counts <- fimo_counts %>%
    rowwise() %>% 
    mutate(
      test_stat = chisq.test(c(down, up),
                             p=c(tbl_directionchange[["down"]] / setsize,
                                 tbl_directionchange[["up"]] / setsize) 
      )$statistic,
      p_val = chisq.test(c(down, up),
                         p=c(tbl_directionchange[["down"]] / setsize,
                             tbl_directionchange[["up"]] / setsize) 
      )$p.value
    )
  fimo_counts$p_adj <-p.adjust(fimo_counts$p_val, method="fdr") 
  # available methods: c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
  
  # plot fimo results
  #------------------------------------------------
  colnames(fimo_counts) <- c("motif_name","counts_down","ns","counts_up","test_stat","p_val","p_adj")
  
  fimo_counts <- fimo_counts %>%
    mutate(
      categories = case_when(grepl("NR3C", motif_name, ignore.case = TRUE) | motif_name=="Ar" | motif_name=="Nr2F6" ~ "NR",
                             grepl("NFKB", motif_name, ignore.case = TRUE) | motif_name=="REL" | motif_name=="RELA" | motif_name=="RELB" ~ "NFKB",
                             grepl("CREB", motif_name, ignore.case = TRUE) | motif_name=="Atf1"| motif_name=="CREM" ~ "CREB",
                             grepl("POU", motif_name, ignore.case = TRUE)  ~ "POU",
                             grepl("KLF", motif_name, ignore.case = TRUE)  ~ "KLF",
                             grepl("JDP", motif_name, ignore.case = TRUE) ~ "AP-1",
                             grepl("JUN", motif_name, ignore.case = TRUE) ~ "AP-1",
                             grepl("IRF", motif_name, ignore.case = TRUE) ~ "IRF",
                             grepl("Stat", motif_name, ignore.case = TRUE) ~ "STAT",
                             TRUE ~ "Other")
    )
  
  # compute fit and plot regression line
  fit <- lm(fimo_counts$counts_up ~ fimo_counts$counts_down, data = fimo_counts)
  gg <- ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point(colour="grey") +
    stat_smooth(method = "lm", col = "blue") +
    
    ggrepel::geom_label_repel(data= fimo_counts %>% dplyr::arrange(p_adj) %>% head(20) %>% dplyr::filter(counts_down<counts_up) ,
                              aes(x=counts_down, y=counts_up, label=motif_name, colour=categories), size=2, 
                              min.segment.length = 0,
                              position = ggrepel::position_nudge_repel(x=-(seqwidth/2), y=(seqwidth/2)),
                              label.padding = 0.1, box.padding = 0.1, show.legend=FALSE)+

    ggrepel::geom_label_repel(data= fimo_counts %>% dplyr::arrange(p_adj) %>% head(20) %>% dplyr::filter(counts_down>counts_up) ,
                              aes(x=counts_down, y=counts_up, label=motif_name, colour=categories), size=2,
                              min.segment.length = 0, force_pill=0.5, force = 5, max.overlaps = 50,
                              position = ggrepel::position_nudge_repel(x=(seqwidth/2), y=-(seqwidth/2)),
                              label.padding = 0.1, box.padding = 0.1, show.legend=FALSE)+
    
    #highlight the same motifs by coloring the point
    geom_point(data= fimo_counts %>% dplyr::arrange(p_adj) %>% head(20),
               aes(x=counts_down, y=counts_up, colour=categories))+
    scale_colour_manual(name="Motif family",
                        values=c("NR"="firebrick","NFKB"="darkblue","IRF"="seagreen","KLF"="deeppink",
                                 "STAT"="darkolivegreen3","POU"="coral","Other"="black"))+
    
    expand_limits(x=-(seqwidth/2), y=-(seqwidth/2))+ # use the slop to dynamically code this (when windows around summit are smaller, we don't have overplotting issues in the lower left)
    
    labs(title=paste0("Peakregion ",seqwidth," bp."),
         x="#Motifmatches in peaks of downregulated genes",
         y="#Motifmatches in peaks of upregulated genes",
         subtitle = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 2) ) ) +
    theme(legend.position = "bottom")
  results=list()
  results[["fimo_counts"]] <- fimo_counts
  results[["plot"]] <- gg
  return(results)
}

fimo_100_results <- make_fimo_chisquare_plot(summit_flank_100bp,fimo_results, 100)
gg_fimo_100 <- fimo_100_results[["plot"]]
fimo_100_results[["fimo_counts"]] %>% arrange(p_adj) %>% head(n=20)

fimo_1000_results <- make_fimo_chisquare_plot(summit_flank_1000bp,fimo_results, 1000)
gg_fimo_1000 <- fimo_1000_results[["plot"]]

#undebug(make_fimo_chisquare_plot)

#---------------------------------
# --- merge figures into panel
#---------------------------------

gg_r1 <-  ggpubr::ggarrange(gg_distbychange, gg_perm_dist , 
                            labels = c("A","B"), widths = c(1,0.6),
                            ncol = 2, nrow=1)

gg_r2 <-  ggpubr::ggarrange(gg_streme_up, gg_streme_down, 
                            labels = c("D","E"), widths = c(1,1),
                            ncol = 2, nrow=1)

gg_c1 <-  ggpubr::ggarrange(gg_r1, gg_r2, 
                            labels = c(NA, NA),
                            ncol = 1, nrow=2, heights=c(1,1))


full_top <- ggpubr::ggarrange(gg_c1, gg_deeptools, 
                              labels = c(NA, "C"),
                              ncol = 2, nrow=1,
                              widths = c(1.7,1))

gg_bottomrow <-  ggpubr::ggarrange(gg_fimo_100, gg_fimo_1000 , 
                                   labels = c("F","G"), widths = c(1,1),
                                   ncol = 2, nrow=1)


full_panel <- ggpubr::ggarrange(full_top, gg_bottomrow, 
                                labels = c(NA, NA),
                                nrow = 2, ncol=1,
                                heights = c(1.7,1))
full_panel

ggsave(here("./results/current/Figures/Figure_peakgeneannotation.png"), full_panel,
       width=190, height=200, units="mm",
       bg="white")

ggsave(here("./results/current/Figures/Figure_peakgeneannotation.pdf"), full_panel,
       width=190, height=200, units="mm",
       bg="white")
