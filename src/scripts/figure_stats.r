suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))

option_list <- list(
  make_option(c( "--tfactivity"),
              type="character",
              help="TF activity, computed using macrophage specific TSS"),
  make_option(c("--expr"),
              type="character",
              help="Path to file with normalized counts, aggregated per mgi symbol"),
  make_option(c("--difffootprint"),
              type="character",
              help="Path to txt file with differential statistics from footprinting analysis"),
  make_option(c("--memedb_expressed"),
              type="character",
              help="Path to rds file of memedb motifs filtered for those expressed"),
  make_option(c("--heatmap"),
              type="character",
              help="Path to deeptools heatmap of motif of interest"),
  make_option(c("--chipms"),
              type="character",
              help="Path to xlsx file with statistics on GR-ChIPMS analysis")
  )

opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(library(here, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(dplyr, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(universalmotif, warn.conflicts=F, quietly=T))

#set defaults for ggplot2 figures
theme_update(panel.background = element_rect(fill = "transparent", colour = NA),
             plot.background = element_rect(fill = "white", colour = NA),
             legend.background = element_rect(fill = "transparent", colour = NA),
             legend.key = element_rect(fill = "transparent", colour = NA),
             text=element_text(size=10, family = "ArialMT", colour="black"),
             title=element_text(size=8, family="ArialMT", colour="black"),
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

#-------------------------------------------------------------------------------
# Load input files
#-------------------------------------------------------------------------------

meme_db_expressed <- readRDS(opt$memedb_expressed)
expr <- read.delim(opt$expr)
TFA <- read.delim(opt$tfactivity)
footprint <- read.table(opt$difffootprint, header=TRUE)

#-------------------------------------------------------------------------------
# Differential footprints
#-------------------------------------------------------------------------------

footprint <- footprint %>% 
  mutate(motif_name = stringr::str_split_fixed(Motif, "\\.",3)[,3],
         Protection_score_diff = Protection_Score_DexLPS - Protection_Score_LPS,
         categories = case_when(grepl("NR3C", motif_name, ignore.case = TRUE) | motif_name=="Ar" | motif_name=="Nr2F6" ~ "NR",
                                  grepl("CREB", motif_name, ignore.case = TRUE) | motif_name=="Atf1"| motif_name=="CREM" ~ "CREB",
                                  grepl("KLF", motif_name, ignore.case = TRUE)  ~ "KLF",
                                  grepl("ATF", motif_name, ignore.case = TRUE) ~ "AP-1",
                                  grepl("JDP", motif_name, ignore.case = TRUE) ~ "AP-1",
                                  grepl("JUN", motif_name, ignore.case = TRUE) ~ "AP-1",
                                  grepl("IRF", motif_name, ignore.case = TRUE) ~ "IRF",
                                  TRUE ~ "Other"),
         padj = p.adjust(P_values, method = "fdr")
         )

# filter on those that we tested in the GLMs
# even if we do, FDR correction kills almost all significance
#footprint_ftr <- footprint %>% 
#  filter(motifnames %in% meme_db_expressed$altname) %>%
#  mutate( padj = p.adjust(P_values, method = "fdr"))


#should we additionally filter on the protection score or not? If we argue for transient binding, we might not want to
table(footprint %>% filter(Num>100 & P_values < 0.05) %>% dplyr::pull(categories))

gg_difffootprint <- ggplot()+
  geom_point(data = footprint %>% filter(Num>100 ), 
             aes(x=Num, y=-log10(P_values)), colour="grey") +
  geom_point(data = footprint %>% filter(Num>100 & P_values < 0.05), 
             aes(x=Num, y=-log10(P_values), colour=categories))+
  ggrepel::geom_label_repel(data = footprint %>% filter(Num>100 & P_values < 0.05 & categories=="Other") ,
                          aes(x=Num, y=-log10(P_values), label=motif_name, colour=categories), size=2.5, min.segment.length = 0,
                          position = ggrepel::position_nudge_repel(x=8000, y=0.1), 
                          label.padding = 0.1, box.padding = 0.1, show.legend = FALSE)+
  scale_colour_manual(name="Motif family",
                      values=c("NR"="firebrick",
                               "AP-1"="dodgerblue",
                               "CREB"="purple",
                               "IRF"="seagreen",
                               "KLF"="deeppink",
                               "Other"="black"))+
  ylim(c(0,3))

gg_difffootprint
#footprint %>% filter(Num>100 & Protection_Score_DexLPS>1) %>% arrange(P_values) %>% View()
#footprint %>% arrange(desc(abs(Protection_score_diff))) %>% head(20)

# requires poppler-cpp system installation
#magick::image_read_pdf(here("../mac_atacseq/results/current/footprints/DexLPSvLPS_diff_footprint/Lineplots/MA1127.1.FOSB::JUN.pdf"))


gg_placeholder <- ggplot() + 
  theme(plot.margin = margin(t=0, l=0, r=0, b=0, unit = "cm"),
        axis.line = element_blank())

#-------------------------------------------------------------------------------
# Expression levels of STATs
#-------------------------------------------------------------------------------

expr_long <- expr %>% 
  tibble::rownames_to_column(var="mgi_symbol") %>%
  tidyr::pivot_longer(cols=-c(mgi_symbol), # everything except mgi_symbol column
                      names_to="sample") %>%
  tidyr::extract(
    "sample", c("condition","rep","type"), regex = "(.*)([0-9])_([^_]+)$") %>%
  mutate(condition=factor(condition, 
                          levels=c("V", "LPS", "LPS_Dex"),
                          labels=c("Veh", "LPS", "DexLPS")))

genes_OI <- c("Stat1","Stat2","Stat3","Stat4","Stat5a","Stat5b","Stat6","Tcf7","Pou2f1","Rel","Nfkb1","Meis1","Nr3c1")

gg_expr <- 
  ggplot(expr_long %>% filter(mgi_symbol %in% genes_OI) %>%
           mutate(across(mgi_symbol, factor, levels=genes_OI)) )+
  geom_abline(slope=0, intercept = 0, colour="grey")+
  #ggplot(expr_long %>% filter(mgi_symbol %in% c("Stat2","Stat3","Tcf7","Pou2f1","Rel","Nfkb1","Meis1","Nr3c1")))+
  geom_boxplot(aes(x=condition, y=log2(value)))+
  geom_point(aes(x=condition, y=log2(value)))+
  facet_wrap(~mgi_symbol, ncol=7)+
  labs( y="log2 (FPKM)")+
  theme(axis.text.x = element_text(angle=30, hjust=1))

gg_expr

#-------------------------------------------------------------------------------
# TF activities
#-------------------------------------------------------------------------------

TFA_long <- as.data.frame(TFA) %>% 
  tibble::rownames_to_column( "TF") %>% 
  tidyr::pivot_longer(cols=2:(ncol(TFA)+1),
                      names_to="sample",
                      values_to="tfactivity") %>%
  tidyr::extract(
    "sample", c("condition","rep","type"), regex = "(.*)([0-9])_([^_]+)$") %>%
  mutate(condition=factor(condition, 
                          levels=c("V", "LPS", "LPS_Dex"),
                          labels=c("Veh", "LPS", "DexLPS")))

gg_TFA_all <- ggplot(TFA_long) +
  geom_boxplot(aes(x=condition, y=tfactivity))+
  facet_wrap(~TF, scales = "free", ncol=6)

ggsave(here("./results/current/tfactivity/TFA_all.png"), gg_TFA_all,
       width=190, height = 250, units="mm")

gg_tfa <- ggplot(TFA_long %>% filter(TF %in% c("NR3C1","STAT2","STAT3","RELA","JUN","JUNB","JUND","FOS","FOSL"))) +
  geom_boxplot(aes(x=condition, y=tfactivity))+
  geom_point(aes(x=condition, y=tfactivity))+
  labs( y = "TF activity")+
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  facet_wrap(~TF, scales = "free")

#-------------------------------------------------------------------------------
# ChIP-MS data
#-------------------------------------------------------------------------------

chipMS <- openxlsx::read.xlsx( opt$chipms )
# save original column names to know what is what
orig_cnames <- colnames(chipMS)
# make them valid for R
colnames(chipMS) <- make.names(colnames(chipMS))
#orig_cnames
#colnames(chipMS)

genes_OI <- c("Rela","Rel","Junb","Nfkb1","Nr3c1","Meis1","Stat1","Stat2","Stat3","Stat4","Stat5b;Stat5a","Stat6")

gg_chipMS <- ggplot(chipMS, aes(x=Test.statistic ,y=X.Log.t.test.p.value ))+
  geom_point(colour="grey")+
  geom_abline(slope=0, intercept = 0, colour="grey")+
  geom_point(data=chipMS %>% filter(X.Log.t.test.p.value > 1.3),
             colour="black")+
  geom_point(data=chipMS %>% filter(Gene.names %in% genes_OI), 
             colour="blue", size=2)+
  ggrepel::geom_label_repel(data=chipMS %>% filter(Gene.names %in% genes_OI), 
                            aes( label=Gene.names ), size=2, nudge_x = 0.5, nudge_y = 0.2)+
  geom_hline(yintercept=1.3, linetype=2)+
  labs(x="Test statistic - wtGR vs wtIgG",
       y="-log(p) - wtGR vs wtIgG")
gg_chipMS

#-------------------------------------------------------------------------------
# deeptools at STAT
#-------------------------------------------------------------------------------

#deeptools <- png::readPNG( opt$heatmap )

#gg_deeptools <- ggplot() + 
#  ggpubr::background_image(deeptools) +
#  # so it doesn't get squished
#  #coord_fixed()+
#  # This ensures that the image leaves some space at the edges
#  theme(plot.margin = margin(t=0, l=0, r=0, b=0, unit = "cm"),
#        axis.line = element_blank())
# --> moved to supplements

#-------------------------------------------------------------------------------
# STAT motifs
#-------------------------------------------------------------------------------

gg_statmotifs <-
  meme_db_expressed %>%
  filter( grepl("STAT", stringr::str_to_upper(altname)) ) %>%
  mutate(name = paste(stringr::str_to_upper(altname), "(", name, ")")) %>%
  to_list() %>% 
  view_motifs()

gg_statmotifs

ggsave(here("./results/current/Figures/Suppl_Figure_statmotifs.pdf"), gg_statmotifs,
       width=180, height=200, units="mm",
       bg="white")

#------------------------------------------------------------------------------------
#------ plot aggregated figure
#------------------------------------------------------------------------------------

#gg_c1 <-  ggpubr::ggarrange(gg_expr, gg_tfa, 
#                            labels = c("A","B"),
#                            ncol = 1, nrow = 2, heights = c(1,1))

left <- ggpubr::ggarrange(gg_expr, gg_difffootprint, 
                          labels = c("A","C"),
                          ncol = 1, nrow=2, heights = c(1,1))

right <- ggpubr::ggarrange(gg_chipMS, gg_placeholder,
                                labels = c("B","D"),
                                ncol = 1, nrow=2, heights = c(2,1))
  
  
full_panel <- ggpubr::ggarrange(left, right, 
                            labels = c(NA,NA),
                            ncol = 2, nrow=1, widths = c(1,0.7))
full_panel


ggsave(here("./results/current/Figures/Figure_stats.png"), full_panel,
       width=190, height=120, units="mm",
       bg="white")

ggsave(here("./results/current/Figures/Figure_stats.pdf"), full_panel,
       width=190, height=120, units="mm",
       bg="white")
