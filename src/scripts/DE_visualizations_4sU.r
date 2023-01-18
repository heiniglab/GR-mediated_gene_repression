suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))

option_list <- list(
  make_option(c("-n", "--norm"),
              type="character",
              help="Path to normalized counts"),
  make_option(c("--contrast"),
              type="character",
              help="Path to annotated tsv file of DeSeq2 contrast"),
  make_option(c("-a", "--annotation"),
              type="character",
              help="Path to metadata with annotation"),
  make_option(c( "--log2fcthresh"),
              type="numeric",
              help="Log2FC threshold used in addition to adj.pval to define significant genes"),
  make_option(c("-o", "--outdir"),
              type="character",
              help="Path of output directory")
    )

opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(library(ComplexHeatmap, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(topGO, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(dplyr, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts=F, quietly=T))

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

dir.create(opt$outdir)
res_DexLPSvLPS_ext <- read.delim( opt$contrast )


#-------------------------------------------------
#---------Check for enrichements
#-------------------------------------------------

# Note: GenTable returns the scores of topGOresult as character which is a problem, since some pvalues are simply returned as "1e-30" and turn into NA when coercing back to numeric
# that's why why made our own custom version of it, where we set the eps cutoff for these small values to FALSE
source("./src/scripts/mytopGOGenTable.R")


maketopGO <- function(ont, all_genes) {
  topgo <- new("topGOdata",
               description = "Test", ontology = ont,
               allGenes = all_genes, nodeSize = 20,
               annot = annFUN.org, mapping="org.Mm.eg.db", ID="ensembl")
  
  # After initializing the topGO object, we can perform the significance tests with a number of differene methods
  resultFisher <- runTest(topgo, algorithm = "classic", statistic = "fisher")
  resultFisher.parentchild <- runTest(topgo, algorithm = "parentchild", statistic = "fisher")
  
  showSigOfNodes(topgo, score(resultFisher.parentchild),
                 firstSigNodes = 10, useInfo ='all')
  
  # initialize empty list for results we want to return later
  results = list()
  #print("Top 50 nodes sorted by Fisher parentchild")
  res_gentable <- mytopGOGenTable(topgo, classicFisher = resultFisher,
                                  parentchildFisher = resultFisher.parentchild,
                                  #classicKS = resultKS, elimKS = resultKS.elim,
                                  orderBy = "parentchildFisher", ranksOf = "classicFisher", topNodes = 50)
  results[["gentable"]] <- res_gentable
  
  # also retrieve and return the gene mapping of those top significant categories
  # get the GO IDs
  top_GOs <- res_gentable %>%
    mutate(parentchildFisher = as.numeric(parentchildFisher)) %>%
    mutate(classicFisher = as.numeric(classicFisher)) %>%
    filter(parentchildFisher<0.05) %>%
    top_n(50, wt=-parentchildFisher) %>%
    pull(GO.ID)
  
  # save the genes within them and their scores to a list
  for (GO in top_GOs){
    df <- data.frame(scores=unlist(scoresInTerm(topgo,GO)),
                     genes=unlist(genesInTerm(topgo,GO)))
    DEgenes <- df %>%
      filter(scores==2)
    
    GOname <- res_gentable %>%
      filter(GO.ID == GO) %>%
      pull(Term)
    GOname <- paste(GO, ":", GOname)
    
    results[[GOname]] <- DEgenes$genes
  }
  
  return(results)
  
}

# We create a binary vector indicating which of the investigates genes came up as significant
all_sig <- as.integer(res_DexLPSvLPS_ext$padj < 0.05)
names(all_sig) <- res_DexLPSvLPS_ext$Row.names
up_sig <- as.integer( (res_DexLPSvLPS_ext$padj < 0.05) & (res_DexLPSvLPS_ext$log2FoldChange > opt$log2fcthresh) )
names(up_sig) <- res_DexLPSvLPS_ext$Row.names
down_sig <- as.integer( (res_DexLPSvLPS_ext$padj < 0.05) & (res_DexLPSvLPS_ext$log2FoldChange < (-opt$log2fcthresh)) )
names(down_sig) <- res_DexLPSvLPS_ext$Row.names

for (cat in c("BP", "MF")){ # iterate over GP categories of interest
  for (set in c("all","down","up")){ #iterate over our genesets
    print(cat)
    print(set)
    if (!file.exists(paste0("results/current/rnaseq_4sU/figures/topGO_enrichment_",cat,"_",set,".rds")) ){
      # open graphics device, that the outgenerated plot will get saved to when calling topGO
      png(paste0(opt$outdir,"/topGO_enrichment_",cat,"_",set,"_network.png"))
      assign(paste0(cat,"_",set), maketopGO(cat, factor(get(paste0(set,"_sig"))) ) )
      dev.off()
      # save the returned object as rds, so it doesn't have to get rerun next time
      saveRDS(get(paste0(cat,"_",set)),
              file=paste0("results/current/rnaseq_4sU/figures/topGO_enrichment_",cat,"_",set,".rds"))
      } else {
        assign(paste0(cat,"_",set), readRDS(file=paste0("results/current/rnaseq_4sU/figures/topGO_enrichment_",cat,"_",set,".rds")) )
        }
  }
}

# We start by looking into the "BP"(=Biological Process) ontology
# pdf(paste0(opt$outdir,"/topGO_enrichment_BP_all_network.pdf"), width=6, height=7, useDingbats = F, pointsize=5)
# BP_all <- maketopGO("BP", factor(all_sig))
# dev.off()
# pdf(paste0(opt$outdir,"/topGO_enrichment_BP_down_network.pdf"), width=6, height=7, useDingbats = F, pointsize=5)
# BP_down <- maketopGO("BP", factor(down_sig))
# dev.off()
# pdf(paste0(opt$outdir,"/topGO_enrichment_BP_up_network.pdf"), width=6, height=7, useDingbats = F, pointsize=5)
# BP_up <- maketopGO("BP", factor(up_sig))
# dev.off()

# And now for the molecular function

# pdf(paste0(opt$outdir,"/topGO_enrichment_MF_all_network.pdf"), width=6, height=7, useDingbats = F, pointsize=5)
# MF_all <- maketopGO("MF", factor(all_sig))
# dev.off()
# pdf(paste0(opt$outdir,"/topGO_enrichment_MF_down_network.pdf"), width=6, height=7, useDingbats = F, pointsize=5)
# MF_down <- maketopGO("MF", factor(down_sig)) 
# dev.off()
# pdf(paste0(opt$outdir,"/topGO_enrichment_MF_up_network.pdf"), width=6, height=7, useDingbats = F, pointsize=5)
# MF_up <- maketopGO("MF", factor(up_sig))
# dev.off()


for (res in c("BP_all","BP_up","BP_down","MF_all","MF_up","MF_down")){
  gg <- get(res)[[1]] %>% 
    mutate(parentchildFisher=as.numeric(parentchildFisher))%>%
    top_n(30, wt=-parentchildFisher) %>% 
    mutate(hitsperc=Significant*100/Annotated) %>% 
    ggplot(aes(x=-log10(parentchildFisher), 
               y=reorder(Term,-parentchildFisher), 
               colour=hitsperc, 
               size=Annotated)) +
    geom_point() +
    scale_colour_gradient(low = "skyblue3", high = "purple2", guide="colourbar",limits=c(0,100))+
    guides(colour = guide_colourbar(barheight = 7))+
    geom_vline(xintercept=-log10(0.05), linetype="dashed")+
    expand_limits(x=0) +
    labs(x="-log10(pvalue)", y="GO term", colour="Hits (%)", size="Termsize")
  ggsave(paste0(opt$outdir,"/topGO_enrichment",res,".png"), gg)
  assign(paste0("gg_enrichment_",res),gg)
}

#-------------------------------------------------
#
#-------------------------------------------------

anno <- read.delim(opt$annotation)

#-------------------------------------------------
#---------library normalization
#-------------------------------------------------

normalized_counts <- read.delim(opt$norm, 
            header = TRUE)

counts_long <- as.data.frame(normalized_counts) %>%
  tidyr::pivot_longer(everything(),
                      names_to="sample",
                      values_to = "counts")

counts_boxplot <- ggplot(data=counts_long)+
  geom_boxplot(aes(x=sample, y=log(counts)))+
  theme(axis.text.x = element_text(angle=90))
counts_boxplot
ggsave(paste0(opt$outdir,"/counts_boxplot.png"), counts_boxplot)
#-------------------------------------------------
#---------PCA
#-------------------------------------------------

# Run PCA on log of normalized counts
#PCA_all = prcomp(log2(normalized_counts+1) %>% t(), 
#                center=TRUE, 
#                scale=FALSE) #center to change mean to 0, scale to change SD to 1


transposed_variables <- as.data.frame(t(log2(normalized_counts+1))) %>%
  mutate(treat=factor(anno$treat, levels=c("Vehicle","LPS","LPS_Dex")) ) %>%
  dplyr::select("treat",everything())

PCA_all <- FactoMineR::PCA(X = transposed_variables, # transposed so that the variables (=genes) are columns
                           scale.unit = TRUE, # whether we scale the variance or not (centering gets done automatically)
                           ncp = 5, # Number of PCs
                           quali.sup = 1, # Position of the qualitative variable
                           graph = FALSE) # No graphic outputs (lib. factoextra)

gg_elbow <- factoextra::fviz_eig(PCA_all, 
                     addlabels = FALSE,
                     ncp=10,
                     main=NULL,
                     ylab="% var",
                     ggtheme = theme(plot.title = element_blank())
)
  

gg_pca1 <- factoextra::fviz_pca_ind(PCA_all,
                         #col.ind="cos2",
                         axes=c(1,2),
                         pointsize=4,
                         habillage = PCA_all$call$X$treat,
                         legend.title = "",
                         title="",
                         #addEllipses = TRUE, #too few points for ellipse
                         label=FALSE,
                         palette = c("#CD8500","#8DA0CB","#66C2A5"), #c("veh" =  "#CD8500", "lps" = "#8DA0CB","lps_dex" = "#66C2A5"))
                         invisible="quali") # in order to disable displaing centroid 

gg_pca2 <- factoextra::fviz_pca_ind(PCA_all,
                         #col.ind="cos2",
                         axes=c(3,4),
                         pointsize=4,
                         habillage = PCA_all$call$X$treat,
                         legend.title = "",
                         #addEllipses = TRUE, #too few points for ellipse
                         label=FALSE,
                         palette = c("#CD8500","#8DA0CB","#66C2A5"), #c("veh" =  "#CD8500", "lps" = "#8DA0CB","lps_dex" = "#66C2A5"))
                         invisible="quali") # in order to disable displaing centroid


gg_pca <- ggpubr::ggarrange(gg_pca1,gg_pca2,
                            common.legend=TRUE)
gg_pca
ggsave(paste0(opt$outdir,"/PCA.png"),gg_pca)

#-------------------------------------------------
#---------Volcano plots
#-------------------------------------------------

res_DexLPSvLPS_ext <- res_DexLPSvLPS_ext %>% 
      mutate(change=case_when(padj<0.05 & log2FoldChange > opt$log2fcthresh ~ "sig up",
                              padj<0.05 & log2FoldChange < (-opt$log2fcthresh) ~ "sig down",
                              TRUE ~ "no change"))
table(res_DexLPSvLPS_ext$change)

gg_volcano <- res_DexLPSvLPS_ext %>%
  ggplot()+
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=change))+
  coord_cartesian(ylim=c(-5,120),xlim=c(-8,8))+
  scale_color_manual(name="",
                     values = c("no change"="#000000","sig down"="#8DA0CB","sig up"="#66C2A5"))+
  ggrepel::geom_text_repel(data=dplyr::filter(res_DexLPSvLPS_ext,-log10(padj) > 28 | abs(log2FoldChange) > 6), 
                           aes(x=log2FoldChange, y=-log10(padj), label=mgi_symbol, colour=change), size=2,
                           min.segment.length = 0, max.overlaps = 100, box.padding = 0.1,
                           show.legend = FALSE)+
  labs(title="DexLPS vs LPS")+
  theme(legend.position = c(0.15, 0.8))
gg_volcano
ggsave(paste0(opt$outdir,"/volcano_plot_DexLPSvLPS.png"), gg_volcano)

#-------------------------------------------------
#---------Heatmap, ALL genes
#-------------------------------------------------

# For the visualizations in heatmaps, we use the normalized counts (counts normalized by DESeq's sizefactors) and add a pseudocount (+1) before taking the log.  
# We then z-scale those log-transformed counts row-wise (for each gene).

LOG.all_zn <- t(apply(log(normalized_counts+1),
                      1, 
                      function(x) (x-mean(x))/sd(x)
                      )
                ) #z-scale normalized values by row (by gene)

#create annotation labels for the heatmap
ha.tmp <- HeatmapAnnotation(df = anno %>% dplyr::select("treat"),
                            col = list(treat = c("Vehicle" =  "#CD8500", 
                                                 "LPS" = "#8DA0CB", 
                                                 "LPS_Dex" = "#66C2A5")),#red,green, blue
                            show_annotation_name = c(bar = FALSE),
                            show_legend = c(treat=TRUE),
                            annotation_legend_param =  list(
                              labels_gp = gpar(fontsize=6),
                              title_gp = gpar(fontsize=6)
                            )
                                       
                            #subset=setNames(RColorBrewer::brewer.pal(6,"Set1") ,levels(met.all$subset))
                            #annotation_legend_param=list( )  #instead of manually assigning each value a color, we can just assign our factor levels the colours of a brewer pallette
                            )


heatmap_allgenes <- Heatmap(
  LOG.all_zn, 
  top_annotation = ha.tmp, 
  column_title="All genes, all samples",
  column_title_gp = gpar(fontsize = 6, fontface = "bold"), 
  name="Row Z-Score", #Title on top of legend
  clustering_distance_rows = "euclidean", 
  clustering_method_rows = "complete", 
  show_row_dend = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_dend_height = unit(1, "cm"),
  row_names_gp = gpar(fontsize = 6),
  heatmap_legend_param = list(
    title_position = "leftcenter-rot",
    legend_direction="vertical",
    labels_gp = gpar(fontsize=6),
    title_gp = gpar(fontsize=6),
    legend_height = unit(2, "cm"),
    at = c(-2, 0, 2), 
    labels = c("-2", "0", "2")
    )
)

pdf(paste0(opt$outdir,"/heatmap_allgenes.pdf"), width=6, height=7, useDingbats = F, pointsize=5)
draw(heatmap_allgenes, heatmap_legend_side="right", merge_legend = TRUE)
dev.off()

gg_heatmap_allgenes <- grid.grabExpr(draw(heatmap_allgenes, heatmap_legend_side="right", merge_legend = TRUE))
gg_heatmap_allgenes

#-------------------------------------------------
#---------Heatmap, SIGgenes
#-------------------------------------------------

#remove those entries in results that have NA as padj and filter for the significant ones
res_DexLPSvLPS_sig <- res_DexLPSvLPS_ext %>%
  filter(!is.na(padj)) %>%
  filter(padj<0.05) %>%
  filter(abs(log2FoldChange) > opt$log2fcthresh)

heatmap_siggenes <- Heatmap(
  LOG.all_zn[res_DexLPSvLPS_sig$Row.names,] , 
  top_annotation = ha.tmp, 
  column_title="DexLPSvsLPS DE genes",
  column_title_gp = gpar(fontsize = 6, fontface = "bold"), 
  name="Row Z-Score", #Title on top of legend
  clustering_distance_rows = "euclidean", 
  clustering_method_rows = "complete", 
  show_row_dend = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_dend_height = unit(1, "cm"),
  row_names_gp = gpar(fontsize = 6),
  heatmap_legend_param = list(
    title_position = "leftcenter-rot",
    legend_direction="vertical",
    labels_gp = gpar(fontsize=6),
    title_gp = gpar(fontsize=6),
    legend_height = unit(2, "cm"),
    at = c(-2, 0, 2), 
    labels = c("-2", "0", "2")
  )
)

pdf(paste0(opt$outdir,"/heatmap_siggenes_DexLPSvsLPS.pdf"), width=6, height=7, useDingbats = F, pointsize=5)
draw(heatmap_siggenes, heatmap_legend_side="bottom")
dev.off()

gg_heatmap_siggenes <- grid.grabExpr(draw(heatmap_siggenes, heatmap_legend_side="right", merge_legend = TRUE))

#---------------------------------------------------------------
#---------Heatmap for DE genes within GO category of interest
#---------------------------------------------------------------
# "GO:0060089 : molecular transducer activity"
# "GO:0001067 : transcription regulatory region nucleic ..."
mapping <- res_DexLPSvLPS_ext %>%
  filter(padj<0.05 & abs(log2FoldChange) > opt$log2fcthresh) %>%
  filter(Row.names %in% MF_all[["GO:0060089 : molecular transducer activity"]]) %>%
  dplyr::select(c(Row.names,mgi_symbol))

heatmap_GO  <- Heatmap(
  LOG.all_zn[mapping$Row.names,], 
  top_annotation = ha.tmp, 
  column_title="Sig genes in GO:0060089 : \n molecular transducer activity",
  column_title_gp = gpar(fontsize = 6, fontface = "bold"), 
  name="Row Z-Score", #Title on top of legend
  clustering_distance_rows = "euclidean", 
  clustering_method_rows = "complete", 
  row_labels = mapping$mgi_symbol,
  show_row_dend = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  column_dend_height = unit(1, "cm"),
  row_names_gp = gpar(fontsize = 5),
  heatmap_legend_param = list(
    title_position = "leftcenter-rot",
    legend_direction="vertical",
    labels_gp = gpar(fontsize=6),
    title_gp = gpar(fontsize=6),
    legend_height = unit(2, "cm"),
    at = c(-2, 0, 2), 
    labels = c("-2", "0", "2")
  )
)
heatmap_GO

pdf(paste0(opt$outdir,"/heatmap_siggenesGO0_DexLPSvsLPS.pdf"), width=6, height=7, useDingbats = F, pointsize=5)
draw(heatmap_GO, heatmap_legend_side="bottom")
dev.off()

gg_heatmap_GO <- grid.grabExpr(draw(heatmap_GO, heatmap_legend_side="right", merge_legend = TRUE))

#------------------------------------------------------------------------------------
#------ include the GO terms 
#------------------------------------------------------------------------------------

GO_network <- png::readPNG("results/current/rnaseq_4sU/figures/topGO_enrichment_MF_all_network.png")

gg_GO_network <- ggplot() + 
  ggpubr::background_image(GO_network) +
  # so it doesn't get squished
  coord_fixed()+
  # This ensures that the image leaves some space at the edges
  theme(plot.margin = margin(t=0, l=0, r=0, b=0, unit = "cm"),
        axis.line = element_blank())



#---------------------------------------------------------------
#---------Pull the figures together
#---------------------------------------------------------------

#lay <- rbind(c(1,2,2),
#             c(3,4,7),
#             c(3,4,7),
#             c(5,5,6),
#             c(5,5,6))

#m1 <- gridExtra::grid.arrange(gg_elbow, gg_pca2, 
#                              gg_heatmap_allgenes, gg_volcano, 
#                              gg_enrichment_MF_all, gg_heatmap_GO0016772,
#                              gg_GO_network,
#                              layout_matrix=lay)


#ggsave("results/current/Figure_2.png",m1, width = 20, height = 30, units = "cm")

gg_r1 <-  ggpubr::ggarrange(gg_elbow, gg_pca1,
                            labels = c("A", NA),
                            ncol = 2, nrow = 1, widths=c(1.4,1))

gg_r2c1 <-  ggpubr::ggarrange(gg_volcano, gg_enrichment_MF_all, gg_GO_network,
                                 labels = c("B", "D","F"),
                                 ncol = 1, nrow = 3, heights = c(1,1,0.5))

gg_r2c2 <-  ggpubr::ggarrange(gg_heatmap_siggenes, gg_heatmap_GO, 
                                 labels = c("C", "E"),
                                 ncol = 1, nrow=2, heights=c(1,1))

gg_r2 <-  ggpubr::ggarrange(gg_r2c1, gg_r2c2, 
                            labels = c(NA, NA),
                            ncol = 2, nrow=1, widths=c(1,1))

full_panel <- ggpubr::ggarrange(gg_r1, gg_r2, 
                                nrow=2, heights=c(0.4,2))

full_panel

ggsave("results/current/Figures/Figure_rnaseq.png", full_panel,
       width=190, height=250, units="mm",
       bg="white")

ggsave("results/current/Figures/Figure_rnaseq.pdf", full_panel,
       width=190, height=250, units="mm",
       bg="white")