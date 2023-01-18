suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))

option_list <- list(
  make_option(c("--model_coefs_joint"),
              type="character",
              help="Path to rds file with model coefficients of joint models"),
  make_option(c("--model_coefs_sep"),
              type="character",
              help="Path to rds file with model coefficients of models tha include enhancers and promoters separately"),
  make_option(c("--auc"),
              type="character",
              help="Path to rds file with auc results of all models"),
  make_option(c("--motifcounts_summitregion"),
              type="character",
              help="Path to rds file with motifcounts within summitregions"),
  make_option(c("--raw_counts"),
              type="character",
              help="Path to rds file of raw counts within the prox based model")
  )

opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(library(ggplot2, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(dplyr, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(here, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(ggcorrplot, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(RColorBrewer, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(circlize, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(ggpubr, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(ComplexHeatmap, warn.conflicts=F, quietly=T))

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

#----------------------------------------------
#------ load data
#----------------------------------------------

model_coefs_joint <- readRDS ( opt$model_coefs_joint )
model_coefs_sep <- readRDS ( opt$model_coefs_sep ) 
AUC_metrics <- readRDS ( opt$auc ) 

AUC_metrics <- AUC_metrics %>% 
  mutate_at(1:6,as.factor) %>%
  mutate(motifdata = dplyr::recode_factor(motifdata, 
                                          motifcounts_abcregion  = "abcregion",
                                          motifcounts_summitregion = "summitregion")) # make factor names shorter

AUC_metrics <- AUC_metrics %>% 
  mutate(names=rownames(.)) # we can't remove the rownames, since we need them for the heatmap annotation later

#---------------------------------------------
#---- show correlation structure
#---------------------------------------------
motifcounts_summitregion <- readRDS( opt$motifcounts_summitregion ) 
#motifdf <- read.table("~/projects/pipeline_ChIP-nexus/results/current/integrate_RNAseq/fimo_featurematrix/merged_matchedpeaks_bygene_30kb_slop50_thresh0.0001_withgenenames.tsv", header=TRUE)
# drop metadata
motifdf <- motifcounts_summitregion[,-c(1:2)]

# determine most variable motifs
motifdf_var <- motifdf %>% dplyr::summarise(across(where(is.numeric), var))
top30perc_var_motif <- sort(motifdf_var, decreasing = T)[1:ceiling(ncol(motifdf_var)*0.3)] %>% names()
# compute and plot their correlation
cor_mat_topvar <- cor(x=motifdf %>% dplyr::select(all_of(top30perc_var_motif)) )
gg_cor <- ggcorrplot(corr = cor_mat_topvar,
                     hc.order = TRUE, # reorders variables according to their correlations
                     outline.col = "white", 
                     colors = c("#6D9EC1", "white", "#E46726"),
                     tl.cex=3)
gg_cor

ggsave(here("./results/current/Figures/Figure_motifcorrelations.png"), gg_cor,
       width=190, height=190, units="mm",
       bg="white")
ggsave(here("./results/current/Figures/Figure_motifcorrelations.pdf"), gg_cor,
       width=190, height=190, units="mm",
       bg="white")
#---------------------------------------------
#---- AUC barplots
#---------------------------------------------
melted_AUC_metrics <- AUC_metrics %>%
  reshape::melt(., measure.vars=c("net_train","net_test"))

# retrieve name of top 25 best performing models (to filter heatmap)
top_models <- melted_AUC_metrics %>% filter(variable=="net_test") %>% arrange(desc(value)) %>% top_n(25) #(with this cutoff the prox based one gets included)
top_models

# overall best performance
max_net_test <- melted_AUC_metrics %>% filter(variable=="net_test") %>% filter(value==max(value)) #%>% summarize(max(value)) %>% as.numeric()
max_net_test

#remove the proximity based model and add it in shape of a reference line instead
prox_reference <- melted_AUC_metrics %>% filter(variable=="net_test") %>% filter(names=="prox") %>% pull(value) %>% as.numeric()
prox_coefs <- model_coefs_joint %>% filter(prox!=0) %>% dplyr::pull(featurename)


# Create the plot
gg_AUC <- ggplot(melted_AUC_metrics %>% filter(variable=="net_test") %>% filter(names!="prox"))+
  geom_bar(aes(x=weight, y=value, fill=condition, alpha=onlymax), 
           colour="black", size=0.5, width=0.9, 
           stat="identity",
           position=position_dodge())+
  geom_hline(aes(yintercept=prox_reference, linetype="reference"), colour="purple")+
  scale_linetype_manual(name="", values = c(2))+
  facet_grid( cols = vars(sepPromEnh,excludepromoters),
              rows= vars(motifdata),
              labeller = labeller(
                motifdata=c("abcregion"="active regions",
                            "summitregion"="GR summitregions"),
                excludepromoters=c("FALSE"="excl.prom: none",
                                   "onlyNONself"="excl.prom: nonself",
                                   "all"="excl.prom: all"),
                sepPromEnh=c("TRUE"="sep prom and enh",
                             "FALSE"= "aggr prom and enh")
              )
  )+
  scale_fill_manual(
    values=c("#339966","#0066CC", "orange"),
    breaks=c("dexlps","lps", "DexLPS-LPS"),
    labels=c("DexLPS","LPS", "\u0394 DexLPS-LPS")
  )+
  scale_alpha_manual(
    values=c(1,0.5),
    breaks=c("FALSE","TRUE"),
    labels=c("1-to-many", "1-to-1")
  )+
  labs(x="Weight features during aggregation",
       y="AUC",
       alpha="mapping"
  )+
  theme(
    legend.position = "bottom")

gg_AUC

#---------------------------------------------
#---- coefficient heatmap
#---------------------------------------------

#AUC_metrics <- AUC_metrics %>% mutate_if(is.character,as.factor)

make_coefficient_heatmap <- function(model_coefs){
  # check what the input was to filter models accordingly later
  if (substitute(model_coefs) == "model_coefs_joint"){
    is_separate=FALSE
  } else if (substitute(model_coefs) == "model_coefs_sep"){
    is_separate=TRUE
  } else {
    stop("The coefficient input matrix doesn't match the expected size")
  }
  
  #assign any 0 coefficients an NA (so they don't get affected by scaling)
  model_coefs[model_coefs==0] <- NA
  # remove the intercept term for plotting
  #model_coefs <- model_coefs %>% filter(featurename!="(Intercept)")
  
  # count how many models have the motif included with a non-zero coef
  model_coefs_ftr <- model_coefs %>%
    mutate( sum = rowSums(across(where(is.numeric)), na.rm = TRUE)) %>%
    mutate( motifhascoef_counts = rowSums(!is.na(.)) ) %>%
    dplyr::select(any_of(
      c("featurename", "sum", "motifhascoef_counts", top_models$names))) %>% # only plot models with best performance
    mutate(motifhascoef_counts_topmodels = rowSums( !is.na( dplyr::select(., 3:ncol(.)) ) )) %>%
    mutate(., across(3:(ncol(.)-1), ~(scale(.) %>% as.vector))) %>% #scale before filtering on spec. factors
    slice_max(motifhascoef_counts_topmodels,n=40)
  
  # remove extra columns before passing it to heatmap
  model_coefs_numeric <- model_coefs_ftr %>% 
    tibble::column_to_rownames("featurename") %>%
    dplyr::select(!c(sum, motifhascoef_counts, motifhascoef_counts_topmodels))
  
  row_ha = rowAnnotation(
    counts = anno_barplot(model_coefs_ftr$motifhascoef_counts_topmodels)
  )
  
  # take annotation from AUC_metrics
  #create annotation labels for the heatmap
  myreds <- brewer.pal(3,"OrRd")
  col_ha <- HeatmapAnnotation(df = AUC_metrics %>% 
                                filter(names %in% names(model_coefs_numeric)) %>%
                                arrange(match(names,names(model_coefs_numeric))) %>%
                                dplyr::select(c(condition,motifdata,excludepromoters, onlymax,weight)),
                              col = list(
                                condition = c("lps" = "#0066CC", "dexlps" = "#339966", "DexLPS-LPS" = "orange"),
                                motifdata=c("summitregion" = "purple", "abcregion" = "chocolate4"),
                                excludepromoters= c("FALSE"=myreds[1], "onlyNONself"=myreds[2], "all"=myreds[3]),
                                onlymax=c("prox"="black", "TRUE"="lightgrey", "FALSE"="darkgrey"),
                                weight=c("FALSE"="lightgoldenrod", "abcscore"="goldenrod")
                              ),
                              AUC = anno_barplot(AUC_metrics %>%
                                                   filter(names %in% names(model_coefs_numeric)) %>%
                                                   arrange(match(names,names(model_coefs_numeric))) %>%
                                                   dplyr::pull(net_test),
                                                 ylim=c(0.5,0.8)),
                              show_annotation_name = c(FALSE,FALSE, FALSE,FALSE,FALSE,TRUE),
                              show_legend = c(group=TRUE),
                              annotation_legend_param =  list(
                                labels_gp = gpar(fontsize=6),
                                title_gp = gpar(fontsize=6)
                              ),
                              simple_anno_size = unit(2,"mm")
  )
  
  mycolorramp <- circlize::colorRamp2(breaks=c(-1.5,0,1.5),
                                      colors=c("blue","white","red"))
  
  heatmap <- Heatmap(
    col=mycolorramp,
    model_coefs_numeric,
    top_annotation = col_ha, 
    right_annotation = row_ha,
    row_split = 4,
    column_title="Model coefficients",
    column_title_gp = gpar(fontsize = 6, fontface = "bold"), 
    column_names_gp = gpar(fontsize = 6), 
    row_names_gp = gpar(fontsize = 6),
    name="z-scaled coefficient", #Title on top of legend
    clustering_distance_rows = "euclidean", 
    clustering_method_rows = "ward.D2", 
    clustering_distance_columns = "euclidean", 
    clustering_method_columns = "ward.D2", 
    show_column_dend = TRUE,
    show_row_dend = TRUE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    row_names_side = "left",
    column_dend_height = unit(0.5, "cm"),
    
    heatmap_legend_param = list(
      title_position = "leftcenter-rot",
      legend_direction="vertical",
      labels_gp = gpar(fontsize=6),
      title_gp = gpar(fontsize=6),
      legend_height = unit(1, "cm"),
      at = c(-1.5, 0, 1.5), 
      labels = c("-1.5", "0", "1.5")
    )
  )
  
  return (
    grid.grabExpr( draw(heatmap, heatmap_legend_side="right", merge_legend = TRUE) )
  )
  
}

gg_model_coefs_joint_heatmap <-  make_coefficient_heatmap(model_coefs_joint)

#---------------------------------------------
#----  coefficients summitranges
#---------------------------------------------

make_summitranges_coefficients_heatmap <- function(model_coefs){
  #assign any 0 coefficients an NA (so they don't get affected by scaling)
  model_coefs[model_coefs==0] <- NA
  
  model_coefs <- model_coefs %>% 
    #filter(featurename!="(Intercept)") %>% # remove the intercept term for plotting
    dplyr::select("featurename","prox" , contains("summitregion")) # select only models based on summitregion
  
  # count how many models have the motif included with a non-zero coef
  model_coefs_ftr <- model_coefs %>%
    mutate( sum = rowSums(across(where(is.numeric)), na.rm = TRUE)) %>% 
    mutate( motifhascoef_counts = rowSums(!is.na(.)) ) %>%
    mutate(., across(3:(ncol(.)-2), ~(scale(.) %>% as.vector))) %>% #scale before filtering on spec. factors 
    filter( motifhascoef_counts >= 15 )
    #filter(prox!=0)
  
  # remove extra columns before passing it to heatmap
  model_coefs_numeric <- model_coefs_ftr %>% 
    tibble::column_to_rownames("featurename") %>%
    dplyr::select(!c( sum, motifhascoef_counts))
  
  row_ha = rowAnnotation(
    counts = anno_barplot(model_coefs_ftr$motifhascoef_counts)
  )
  
  # take annotation from AUC_metrics
  #create annotation labels for the heatmap
  myreds <- brewer.pal(3,"OrRd")
  col_ha <- HeatmapAnnotation(df = AUC_metrics %>% 
                                filter(names %in% names(model_coefs_numeric)) %>%
                                arrange(match(names,names(model_coefs_numeric))) %>%
                                dplyr::select(c(condition,motifdata,excludepromoters, onlymax,weight)),
                              col = list(
                                condition = c("lps" = "#0066CC", "dexlps" = "#339966", "DexLPS-LPS" = "orange"),
                                motifdata=c("summitregion" = "purple", "abcregion" = "chocolate4"),
                                excludepromoters= c("FALSE"=myreds[1], "onlyNONself"=myreds[2], "all"=myreds[3]),
                                onlymax=c("prox"="black", "TRUE"="lightgrey", "FALSE"="darkgrey"),
                                weight=c("FALSE"="lightgoldenrod", "abcscore"="goldenrod")
                              ),
                              AUC = anno_barplot(AUC_metrics %>% 
                                                   filter(names %in% names(model_coefs_numeric)) %>%
                                                   arrange(match(names,names(model_coefs_numeric))) %>%
                                                   dplyr::pull(net_test),
                                                 ylim=c(0.5,0.8)),
                              show_annotation_name = c(FALSE,FALSE, FALSE,FALSE,FALSE,TRUE),
                              show_legend = FALSE,
                              annotation_legend_param =  list(
                                labels_gp = gpar(fontsize=6),
                                title_gp = gpar(fontsize=6)
                              ),
                              simple_anno_size = unit(2,"mm")
  )
  
  mycolorramp <- circlize::colorRamp2(breaks=c(-1.5,0,1.5),
                                      colors=c("blue","white","red"))
  
  heatmap <- Heatmap(
    col=mycolorramp,
    model_coefs_numeric,
    top_annotation = col_ha, 
    right_annotation = row_ha,
    row_split = 4,
    column_title="Model coefficients",
    column_title_gp = gpar(fontsize = 6, fontface = "bold"), 
    column_names_gp = gpar(fontsize = 6), 
    row_names_gp = gpar(fontsize = 6),
    name="z-scaled coefficient", #Title on top of legend
    clustering_distance_rows = "euclidean", 
    clustering_method_rows = "ward.D2", 
    clustering_distance_columns = "euclidean", 
    clustering_method_columns = "ward.D2", 
    show_column_dend = TRUE,
    show_row_dend = TRUE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    row_names_side = "left",
    column_dend_height = unit(0.5, "cm"),
    show_heatmap_legend = FALSE
  )
  
  return (
    grid.grabExpr( draw(heatmap, heatmap_legend_side="right", merge_legend = TRUE) )
  )
  
}

gg_model_summitregioncoefs_joint_heatmap <-  make_summitranges_coefficients_heatmap(model_coefs_joint)

#---------------------------------------------
#---- compare model coefficients of best models
#---------------------------------------------

# check in AUC what is the best model
max_net_test$names #"motifcounts_abcregion_condition_DexLPS-LPS_exclprom_onlyNONself_onlymax_FALSE_sepPromEnh_FALSE_weight_FALSE"

# keep those rows where at least one of the models has a non-zero coefficient
keepRows <- which(rowSums(model_coefs_joint [,c("prox",max_net_test$names)]) != 0 )

model_coefs_long <- model_coefs_joint[keepRows,] %>% 
  dplyr::select(c(featurename,prox,max_net_test$names)) %>%
  dplyr::rename(abc=max_net_test$names)%>%
  tidyr::pivot_longer(.,
                      col=c('prox','abc'),
                      names_to = "model",
                      values_to = "coefficients")
model_coefs_long <- model_coefs_long %>%
  mutate(my_color = case_when(coefficients>0 & model=="prox" ~ "darkred",
                              coefficients<0 & model=="prox" ~ "darkblue",
                              coefficients>0 & model!="prox" ~ "lightred",
                              coefficients<0 & model!="prox" ~ "lightblue"
    
  ))


gg_bestmodels <- ggplot(data=model_coefs_long)+
  geom_bar(aes(x=reorder(featurename,-coefficients), y=coefficients, fill=model),stat="identity", position = position_dodge(width = 0.5),width=0.5)+
  scale_fill_manual(values=c("prox"="purple",
                             "abc" = "brown"),
                      labels=c("prox. based",
                               "ABC based"))+
  labs(x="motifname", y="coefficient")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = c(0.8, 0.8))
gg_bestmodels
#---------------------------------------------
#---- or just display the prox based one
#---------------------------------------------
keepRows <- which(model_coefs_joint [,"prox"] != 0 )
model_coefs_data <- model_coefs_joint[keepRows,] %>% 
  mutate(my_color = case_when(prox>0  ~ "positive",
                              prox<0  ~ "negative"))

gg_proxmodel <- ggplot(data= model_coefs_data)+
  geom_bar(aes(x=reorder(featurename,-prox), y=prox, fill=my_color), stat="identity", 
           position = position_dodge(width = 0.5),width=0.5)+
  scale_fill_manual(values=c("positive"="red",
                             "negative" = "blue"))+
  labs(x="motifname", y="coefficient", fill="coefficient sign")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position = c(0.8,0.8),
        panel.grid.major = element_blank())
gg_proxmodel

#---------------------------------------------
#---- with red and blue, paired dark and light for the 2 models
#---------------------------------------------

#my_pal <- RColorBrewer::brewer.pal(6,"Paired")
#gg_bestmodels <- ggplot(data=model_coefs_long)+
#  geom_bar(aes(x=featurename, y=coefficients, fill=my_color, colour=model),stat="identity", position = position_dodge(width = 0.5),width=0.5)+
#  scale_fill_manual(values=c("darkred"=my_pal[6],
#                             "darkblue"=my_pal[2],
#                             "lightred"=my_pal[5],
#                             "lightblue"=my_pal[1]),
#                    labels=c("only peakregions\n- proxbased - up", 
#                             "only peakregions\n- proxbased - down", 
#                             "all enhancers\n- abcbased - up",
#                             "all enhancers\n- abcbased - down"
#                             ))+
#  scale_colour_manual(values=c("prox"="red",
#                               "motifcounts_peak_condition_dexlps_exclprom_all_onlymax_FALSE_sepPromEnh_FALSE_weight_abcscore" = "blue"),
#                      labels=c("prox. based",
#                               "ABC based"))+
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#  guides(fill="none")


#---------------------------------------------
#---- load raw counts ggplot
#---------------------------------------------

gg_rawcounts <- readRDS(opt$raw_counts)

#---------------------------------------------
#---- put figure together
#---------------------------------------------



gg_coef_heatmaps <-  ggpubr::ggarrange(gg_model_coefs_joint_heatmap, gg_model_summitregioncoefs_joint_heatmap,
                                       labels = c("B", "C"),
                                       ncol = 2, nrow = 1, widths=c(1.2,1)
                                       )
gg_coef_heatmaps

#gg_bottomrow <-  ggpubr::ggarrange(gg_proxmodel, gg_rawcounts,
#                                       labels = c("D", "E"),
#                                       ncol = 2, nrow = 1, widths=c(2,1)
#)
#gg_bottomrow
#full_panel <- ggpubr::ggarrange(gg_AUC, gg_coef_heatmaps, gg_bottomrow,

full_panel <- ggpubr::ggarrange(gg_AUC, gg_coef_heatmaps, gg_bestmodels,
                                labels = c("A", NA, "D"),
                                nrow=3, heights=c(0.5,1,0.3))

full_panel

ggsave(here("./results/current/Figures/Figure_GLMs.png"), full_panel,
       width=190, height=240, units="mm",
       bg="white")
ggsave(here("./results/current/Figures/Figure_GLMs.pdf"), full_panel,
       width=190, height=240, units="mm",
       bg="white")
#---------------------------------------------
#---- save models with enhancers and promoters as separate features
#---------------------------------------------

# possible supplemental:
gg_model_coefs_sep_heatmap <-  make_coefficient_heatmap(model_coefs_sep)
gg_model_coefs_sep_heatmap


gg <-  ggpubr::ggarrange(gg_model_coefs_sep_heatmap, 
                         labels = c(NA),
                         ncol = 1, nrow = 1, widths=c(1))
gg
ggsave(here("./results/current/Figures/Figure_GLMs_coefssep.png"), gg,
       width=190, height=190, units="mm",
       bg="white")
ggsave(here("./results/current/Figures/Figure_GLMs_coefssep.pdf"), gg,
       width=190, height=190, units="mm",
       bg="white")