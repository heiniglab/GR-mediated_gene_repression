suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))

option_list <- list(
  make_option(c("--summitAnno"),
              type="character",
              help="Path to dataframe object of GR summits annotated to closest gene"),
  make_option(c("--chipseq_peaks"),
              type="character",
              help="Path to IDR ChIPseq peaks"),
  make_option(c("--chipseq_summits"),
              type="character",
              help="Path to summit file of IDR peaks"),
  make_option(c("--sm_summitranges"),
              type="character",
              help="Path to rds file of genomation score matrix around peak summits"),
  make_option(c("--nr3c1fullsitematches"),
              type="character",
              help="Path to homer hits of nr3c1 fullsite motif"),
  make_option(c("--nr3c1halfsitematches"),
              type="character",
              help="Path to homer hits of nr3c1 halfsite motif"),
  make_option(c( "--streme"),
              type="character",
              help="Path to XML streme output file of enrichment around peak summits"))

opt <- parse_args(OptionParser(option_list=option_list))

# set output for logfile to retrieve stats for plot later
sink(file="results/current/figure_chipseq.out")
    
suppressPackageStartupMessages(library(ggplot2, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(grid, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(dplyr, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(plyranges, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(memes, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(universalmotif, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(here, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(rtracklayer, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(ChIPseeker, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(genomation, warn.conflicts=F, quietly=T))

#set defaults for ggplot2 figures
theme_update(panel.background = element_rect(fill = "transparent", colour = NA),
             plot.background = element_rect(fill = "transparent", colour = NA),
             legend.background = element_rect(fill = "transparent", colour = NA),
             legend.key = element_rect(fill = "transparent", colour = NA),
             text=element_text(size=12, family = "ArialMT", colour="black"),
             title=element_text(size=16, family="ArialMT", colour="black"),
             panel.grid.major = element_line(colour="grey", size=0.2),
             panel.grid.minor = element_blank(),
             axis.text = element_text(size=12, family="ArialMT", colour="black"),
             axis.line = element_line(colour="black"),
             axis.ticks = element_line(colour="black"),
             legend.key.size = unit(12, 'points'), #change legend key size
             legend.key.height = unit(12, 'points'), #change legend key height
             legend.key.width = unit(12, 'points'), #change legend key width
             legend.text = element_text(size=8, family="ArialMT", colour="black")
)

#----------------------------------------------
#------ load input data
#----------------------------------------------
chipseq <- rtracklayer::import.bed(opt$chipseq_peaks)
chipseq <- GenomeInfoDb::keepStandardChromosomes(chipseq, pruning.mode = "tidy")
names(chipseq) <- c(1:length(chipseq))
print(paste("We find a total of ", length(chipseq),"ChIPseq peaks"))

chipseq_summits <- read.table(opt$chipseq_summits)
chipseq_summits <- GRanges(seqnames = chipseq_summits[,c("V1")],
                           ranges = IRanges(start=chipseq_summits[,c("V2")],
                                            end=chipseq_summits[,c("V3")]-1)) # to make up for 0 vs 1 encoding
chipseq_summits$id <- c(1:length(chipseq_summits))

chipseq_summitranges <- chipseq_summits %>% 
  plyranges::anchor_center() %>% 
  plyranges::mutate(width = 1000) 

sm_summitranges <- readRDS(opt$sm_summitranges)

#----------------------------------------------
#------ peak width
#----------------------------------------------
xbreaks <- c(200,500,1000,3000)
gg_peakwidth <- 
  ggplot(as.data.frame(chipseq))+
  geom_histogram(aes(log10(width)), bins=50,
               colour="black")+
  scale_x_continuous("width (bp)", breaks=log10(xbreaks), labels=xbreaks, limits=log10(c(200,4000)) )
gg_peakwidth

print(paste("The peak width has a mean of", format(mean(width(chipseq)), digits = 6), 
            "and a median of", format(median(width(chipseq)),digits = 6)))

#----------------------------------------------
#------ peak wrt reference
#----------------------------------------------
summitAnno <- readRDS(opt$summitAnno)

## 1. where do the peaks lie
#-------------------------------
annostat <- as.data.frame(summitAnno@annoStat)
gg_annopie <- ggplot(annostat, aes(x="", y=Frequency, fill=Feature))+
  geom_bar(stat="identity", width=1, color="white") +
  geom_text(aes(label = paste(round(Frequency, 1), "%"), x = 1.5),
            position = position_stack(vjust = 0.5), size=2) +
  coord_polar("y", start=0)+
  #scale_fill_brewer(palette="Set3")+
  scale_fill_manual(values= c(rev(RColorBrewer::brewer.pal(9,"YlGnBu")),"white", "darkseagreen"))+
  theme_void()+
  theme(
    text=element_text(size=6, family = "ArialMT", colour="black"),
    legend.key.height = unit(6, 'points'), #change legend key height
    legend.key.width = unit(6, 'points'), #change legend key width
    legend.text = element_text(size=6, family="ArialMT", colour="black")
  )
gg_annopie

print(paste("Frequency of peaks within promoters <3kb:", format(sum(annostat[1:3,"Frequency"]),digits=4)))
print(paste("Frequency of peaks in introns:", format(sum(annostat[8:9,"Frequency"]),digits=4)))
print(paste("Frequency of peaks in distal intergenic regions:", format(sum(annostat[11,"Frequency"]),digits=4)))


## 2. how far are they from closest TSS
#-------------------------------
xbreaks <- c(-100, -75, -50, -25,  0, 25, 50, 75, 100)
gg_distexpr <- 
  as.data.frame(summitAnno) %>%
  ggplot( aes(x=distanceToTSS)) +
  geom_histogram(binwidth = 3000) +
  coord_cartesian(
    #  ylim=c(0,2000),
    xlim=c(-100*10^3, 100*10^3)
  )+
  geom_segment(aes(x = 30*1000, y = 0, xend = 30*1000, yend = 3000), colour="black", linetype=2)+
  geom_segment(aes(x = -30*1000, y = 0, xend = -30*1000, yend = 3000), colour="black", linetype=2)+
  scale_x_continuous(
    breaks = xbreaks*1000,
    labels = paste(xbreaks, "kb")
  )+
  labs(title="",
       y="Counts")+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
gg_distexpr

print(paste("Out of all our peaks, a total of",
as.data.frame(summitAnno) %>% filter(abs(distanceToTSS) < 30000) %>% nrow(),
"where within 30000kb of a TSS"))

print(paste("This corresponds to",
format(as.data.frame(summitAnno) %>% filter(abs(distanceToTSS) < 30000) %>% nrow() / length(chipseq) *100, digits=4), "%"))

print(paste("They annotated to ",
as.data.frame(summitAnno) %>% filter(abs(distanceToTSS) < 30000) %>% pull(SYMBOL) %>% unique() %>% length(),
"unique gene symbols"))

#----------------------------------------------
#------ load motif query and match
#----------------------------------------------


load_sitematches <- function(bedpath){
  sitematches <- rtracklayer::import.bed(bedpath)
  # remove weird repetition of seqname
  seqlevels(sitematches) <- gsub(" .*$","",seqlevels(sitematches))
  sitematches <- GenomeInfoDb::keepStandardChromosomes(sitematches, pruning.mode = "tidy")
  names(sitematches) <- c(1:length(sitematches))
  return(sitematches)
}

nr3c1fullsitematches <- load_sitematches(opt$nr3c1fullsitematches)
nr3c1halfsitematches <- load_sitematches(opt$nr3c1halfsitematches)

print(paste("Genome-wide,  we found", length(nr3c1fullsitematches), "matches for the fullsite and",
            length(nr3c1halfsitematches), "matches for the halfsite"))

#----------------------------------------------
#------ # counts for overlap with nr3c1 motif
#----------------------------------------------

# crashes if including halfsites
#overlap_wfullsite <- ChIPpeakAnno::findOverlapsOfPeaks(chipseq,
#                                                       nr3c1fullsitematches,
#                                                       nr3c1halfsitematches,
#                                                       minoverlap=1)
#gg_venn <- grid.grabExpr(
#  ChIPpeakAnno::makeVennDiagram(overlap_wfullsite,
#                                fill=c("#669933", "#ff9900", "#c01311"), # circle fill color (green, orange)
#                                col=c("#669933", "#ff9900", "#c01311"), #circle border color
#                                cat.col=c("#669933", "#ff9900", "#c01311"),
#                                method = NULL,
#                                cex = 0.6,
#                                connectedPeaks = "keepFirstListConsistent"),
#  vp = viewport(w = .6, h = 1.0)
#)


chipseq_nr3c1fullsitematches_counts <- plyranges::count_overlaps(chipseq, nr3c1fullsitematches) %>% 
  as.data.frame() %>% 
  magrittr::set_colnames(c("sitecounts")) %>% 
  count(sitecounts) %>%
  mutate(freq = n / sum(n),
         site = "full")

print(paste("The number of peaks which have at least 1 fullsite match:", 
            sum(chipseq_nr3c1fullsitematches_counts$n) - chipseq_nr3c1fullsitematches_counts$n[1] ))
print(paste("This corresponds to",
            sum(chipseq_nr3c1fullsitematches_counts$freq) - chipseq_nr3c1fullsitematches_counts$freq[1], "%" ))

chipseq_nr3c1halfsitematches_counts <- plyranges::count_overlaps(chipseq, nr3c1halfsitematches)%>% 
  as.data.frame() %>% 
  magrittr::set_colnames(c("sitecounts")) %>% 
  count(sitecounts) %>%
  mutate(freq = n / sum(n),
         site="half")

print(paste("The number of peaks which have at least 1 halfsite match:", 
            sum(chipseq_nr3c1halfsitematches_counts$n) - chipseq_nr3c1halfsitematches_counts$n[1] ))
print(paste("This corresponds to",
            sum(chipseq_nr3c1halfsitematches_counts$freq) - chipseq_nr3c1halfsitematches_counts$freq[1], "%" ))

rm(nr3c1fullsitematches)
rm(nr3c1halfsitematches)

chipseq_nr3c1sitematches_counts <- rbind(chipseq_nr3c1fullsitematches_counts,chipseq_nr3c1halfsitematches_counts)
chipseq_nr3c1sitematches_counts <- chipseq_nr3c1sitematches_counts %>%
  mutate(sitecounts_fac = case_when(sitecounts<=15 ~ as.character(sitecounts),
                                    sitecounts>15 ~ ">15"), # aggregate it towards to top end
         sitecounts_fac = factor(as.character(sitecounts_fac)) )
# set order of factor levels
my_levels <- levels(chipseq_nr3c1sitematches_counts$sitecounts_fac)[gtools::mixedorder(levels(chipseq_nr3c1sitematches_counts$sitecounts_fac))]         
chipseq_nr3c1sitematches_counts <- chipseq_nr3c1sitematches_counts %>%
  mutate(sitecounts_fac = factor(sitecounts_fac, levels=my_levels)) %>%
  group_by(site, sitecounts_fac) %>%
  summarize(freq=sum(freq),
            total=sum(n))

gg_barplotmotifhits <- ggplot(data=chipseq_nr3c1sitematches_counts, aes(x=sitecounts_fac, y= freq*100, group=site))+
  geom_bar(aes(colour=site, fill=site), stat="identity",
           alpha=0.7, position = position_dodge2(width=0.4, padding=0.1, preserve = "single") )+
  scale_fill_manual("",
                    labels = c("NR3C1 fullsite", "NR3C1 halfsite"),
                    breaks = c("full", "half"),
                    values = c("black", "darkgrey")) +
  scale_colour_manual("",
                      labels = c("NR3C1 fullsite", "NR3C1 halfsite"),
                      breaks = c("full", "half"),
                      values = c("black", "darkgrey"))+
  labs(x="# of motifmatches",
       y="% of ChIPseq peaks")+
  theme(legend.position = c(0.7, 0.8),
        axis.text.x = element_text(angle=45, vjust=1, hjust=1))
gg_barplotmotifhits


# can we use the scorematrix function from genomation to get a distribution of reads wrt handpicked motifs

#------------------------------------------------------------------------------------
#------ memes
#------------------------------------------------------------------------------------
# conda activate py_3
# perlbrew use perl-5.34.0
# nohup Rscript memes_runanalyses.R & (from within the script directory)

#------ STREME
#-----------------
streme_results <- memes::importStremeXML( opt$streme )
streme_results <- streme_results %>% 
  mutate(name = paste(consensus, pval))

#DT::datatable(streme_results %>% dplyr::relocate("eval",.after="consensus"),
#              extensions = c('Buttons'), 
#              width=1080,
#              options = list(dom = 'Bfrtip', buttons = c('csv', 'excel'), scrollX = TRUE), 
#              caption = htmltools::tags$caption(style = 'caption-side: bottom; text-align: center;', htmltools::em('Motif enrichment - Discriminative analysis'))
#)

gg_streme <-
  streme_results[1:5,] %>% 
  to_list() %>% 
  view_motifs(names.pos = "top", 
              tryRC = FALSE # we don't care about maximizing based on alignment score, just wanna display the discovered motifs
              )
gg_streme

# why are the lists with pos distr not the same length? Let's pad them
padna <- function(myvector, outlength){
  # check length of vector and how was it's off from desired length
  tofill <- outlength - length(myvector)
  if (tofill<0) {
    stop("Desired length is longer than current length")
  } else if(tofill > 1) {
    # check if it's even
    if((tofill %% 2) == 0) {
      outvector = c(rep(NA,tofill/2), myvector,rep(NA,tofill/2))
    } else {
      outvector = c(rep(NA,floor(tofill/2)), myvector,rep(NA,ceiling(tofill/2)))
    }
    # if uneven, put -1 before and +1 after
    return(outvector)
  }
}

pos_distr <- streme_results$site_distr %>%
  stringr::str_trim(side = c("both")) %>%
  stringr::str_split(pattern=" ")
  
pos_distr <- lapply(pos_distr, function(x) padna(x,outlength=101))

pos_distr <- as.data.frame(do.call(rbind, pos_distr))
colnames(pos_distr) <- 1:ncol(pos_distr)

pos_distr_top_long <- pos_distr %>%
  mutate(altname = streme_results$altname) %>%
  filter(altname %in% paste0("STREME-",c(1,2,3,4,5)) ) %>%
  relocate(altname) %>%
  tidyr::pivot_longer(cols=1:ncol(pos_distr)+1, names_to="position", values_to = "signal") %>%
  mutate(position=as.numeric(position),
         signal=as.numeric(signal)) 

gg_pos_distr <- ggplot(pos_distr_top_long,aes(x=position, y=signal))+
  geom_point(size=0.5)+
  geom_smooth()+
  facet_wrap(~altname, ncol=1, scales="free_y")+
  labs(x="potision",y="frequency")+
  scale_x_continuous(breaks=c(1,51,101), labels=c(-50,0,50))+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.ticks.y  = element_blank(),
    axis.text.y  = element_blank()
  )

gg_pos_distr

#------------------------------------------------------------------------------------
#------ read distribution around peak summit ?
#------------------------------------------------------------------------------------

genomation_profiledata <- genomation::plotMeta(sm_summitranges)

# make sure the xcoords match the values specified in the summitranges (or the xaxis labels will be wrong)
gg_chipseq_genomationprofileplot <-
  ggplot()+
  geom_line(aes(x=seq(1:length(genomation_profiledata)),y=genomation_profiledata))+
  scale_x_continuous("bases", breaks=c(0,250,500,750,1000), labels=c(-500,-250,0,250,500))+
  labs(x="bases",y="average score")+
  theme(plot.margin = margin(1,2,0,0, "cm"))

sm_scaled = genomation::scaleScoreMatrix(sm_summitranges)
gg_chipseq_genomationheatmap <- grid.grabExpr(
  genomation::heatMatrix(sm_scaled, xcoords = c(-500, 500)) 
)

#------------------------------------------------------------------------------------
#------ plot aggregated figure
#------------------------------------------------------------------------------------

gg_c1 <-  ggpubr::ggarrange(gg_barplotmotifhits,gg_peakwidth, 
                            labels = c("A","B"),
                            ncol = 1, nrow = 2, heights=c(1,0.8))

gg_c2 <-  ggpubr::ggarrange(gg_chipseq_genomationprofileplot, gg_chipseq_genomationheatmap, 
                            labels = c("C",NA),
                            ncol = 1, nrow = 2, heights=c(1,2))

gg_c3_r1 <-  ggpubr::ggarrange(gg_streme, gg_pos_distr,
                            labels = c("D",NA),
                            ncol = 2, nrow = 1, widths = c(1,0.5)) 

gg_c3_r2 <-  ggpubr::ggarrange(gg_annopie, gg_distexpr,
                            labels = c("E", "F"),
                            ncol = 2, nrow = 1, widths = c(1,1)) 

gg_c3 <-  ggpubr::ggarrange(gg_c3_r1, gg_c3_r2,
                               labels = c(NA,NA),
                               ncol = 1, nrow = 2, heights=c(1,1)) 

full_panel <-  ggpubr::ggarrange(gg_c1, gg_c2, gg_c3,
                                 labels = NA,
                                 ncol = 3, nrow=1, widths = c(1,1,2))
full_panel
# usually 190(width) by 100(height)
# scale it up, so motifs are displayed correctly
ggsave(here("./results/current/Figures/Figure_chipseq.png"), full_panel,
       width=380, height=200, units="mm",
       bg="white")
ggsave(here("./results/current/Figures/Figure_chipseq.pdf"), full_panel,
       width=380, height=200, units="mm",
       bg="white")
sink()
