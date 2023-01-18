suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))

option_list <- list(
  make_option(c( "--contrast_DexLPSvLPS"),
              type="character",
              help="Path to annotated tsv file of DeSeq2 contrast of DexLPS vs LPS"),
  make_option(c("--assignment_summit_prox"),
              type="character",
              help="Path to rds file of proximity based assignment of peak summits to genes"),
  make_option(c("--assignment_summits_abcregion_dexlps"),
              type="character",
              help="Path to rds file of assignment of peak summits within abcregions to genes (in DexLPS condition)"),
  make_option(c("--assignment_summits_abcregion_lps"),
              type="character",
              help="Path to rds file of assignment of peak summits within abcregions to genes (in LPS condition)"),
  make_option(c("--assignment_abcregion_dexlps"),
              type="character",
              help="Path to rds file of assignment of abcregions to genes (in DexLPS condition)"),
  make_option(c( "--assignment_abcregion_lps"),
              type="character",
              help="Path to rds file of assignment of abcregions to genes (in LPS condition)"),
  make_option(c("--motifcounts_summitregion"),
              type="character",
              help="Path to rds file of fimo motifcounts within summitregions"),
  make_option(c("--motifcounts_abcregion_dexlps"),
              type="character",
              help="Path to rds file of fimo motifcounts within ABC regions (in DexLPS condition)"),
  make_option(c("--motifcounts_abcregion_lps"),
              type="character",
              help="Path to rds file of fimo motifcounts within ABC regions (in LPS condition)"),
  make_option(c( "--outdir"),
              type="character",
              help="Path to output directory")
  )

opt <- parse_args(OptionParser(option_list=option_list))

dir.create(opt$outdir)

#change default for stringAsFactors
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages(library(here, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(plyranges, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(org.Mm.eg.db, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(DESeq2, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(ComplexHeatmap, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(dplyr, warn.conflicts=F, quietly=T))

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
             legend.text = element_text(size=6, family="ArialMT", colour="black"))

set.seed(12345)

#-------------------------------
## read in  data
#-------------------------------

contrast_DexLPSvLPS <- read.delim(opt$contrast_DexLPSvLPS)

for (optname in names(opt)[2:9]){ #except for the outdir
  print(paste0("Loading ", optname))
  assign(optname, readRDS( opt[[optname]] ))
}

#--------------------------------------------
## ---- function definitions
#--------------------------------------------

coerce_coef2df <- function(model_coef){
  model_coef <- as.matrix(model_coef)
  model_coef <- as.data.frame(model_coef)
  model_coef$names <- rownames(model_coef)
  colnames(model_coef) <- c("estimates", "names")
  return(model_coef)
}

#------------function will------------:
# * merge the motifdata with gene assignments and 
# * then merge the gene expression changes,
# * filter for DE genes and aggregate per gene
#--------------------------------------
merge_motifdata_with_assignments <- function( motifcounts, assignments, contrast, 
                                              maxonly=FALSE, excludepromoters=FALSE, weightby=FALSE, sepPromEnh=FALSE){
  
  # check if we should use all abc assignments, or only the max one of each peakID
  if(maxonly==TRUE){
    assignments <- assignments %>% group_by(name) %>% filter(abcscore==max(abcscore)) %>% distinct()
  } else{
    assignments <- assignments
  }
  
  if(excludepromoters=="all"){
    assignments <- assignments %>% filter(!class=="promoter")
  } else if (excludepromoters=="onlyNONself"){
    assignments <- assignments %>% filter(!c(class=="promoter" & isSelfPromoter=="False"))
  } else{
    assignments <- assignments
  }
  
  motifdf <- merge(motifcounts, assignments, 
                   by.x="name", by.y="name")
  
  motifdf <- motifdf %>% relocate(c(anno))
  
  # merge the expression change
  # use mgi_symbol for prox based, otherwise ensemblID
  # We DONT set all.x=TRUE because we don't care about predicting gene that aren't even expressed or that we don't have a clear label for
  if(maxonly=="prox"){
    motifdf <- merge(motifdf, contrast, 
                     by.x="anno", by.y="mgi_symbol")
  } else {
    motifdf <- merge(motifdf, contrast, 
                     by.x="anno", by.y="Row.names")
  }
 
  #recode the logFC and padj into a label (optionally through command line arguments)
  motifdf <- motifdf %>% 
    mutate(label=case_when(log2FoldChange>0.58 & padj < 0.05 ~ "up",
                           log2FoldChange<(-0.58) & padj < 0.05 ~ "down",
                           TRUE ~ "no_change")) %>%
    filter(label!="no_change") %>%
    mutate(label=factor(label,
                        levels=c("down","up"),
                        labels=c(0,1))) %>%
    relocate(label)
  
  # chr 2, 3 and 4 (20%) were used as the tuning set for hyperparameter tuning. 
  # Regions from chromosomes 1, 8 and 9 (20%) were used as the test set for performance evaluation 
  # The remaining regions were used for model training.
  
  #------aggregate over gene SYMBOL
  if ("abcscore" %in% colnames(motifdf)) { # loop for ABC based assignments
    
    if (weightby=="abcscore"){
      unselect_col <- "abcnumerator"
    } else if(weightby=="abcnumerator") {
      unselect_col <- "abcscore"
    } else{
      unselect_col <- c("abcscore","abcnumerator")
    }
    
    if (sepPromEnh==TRUE){
      motifdf_aggr <-
        motifdf %>% 
        dplyr::select(!c(unselect_col,"name","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","gene_biotype","mgi_symbol")) %>%
        { if(weightby!=FALSE) mutate(.,across(where(is.numeric), ~ (.x * get(weightby)))) else .} %>% # weight features by score
        dplyr::select(!any_of(as.character(weightby))) %>% # then we can drop the score since it will be nonsensical after the aggregation anyways
        group_by(label,seqnames,anno,class) %>%
        dplyr::summarise(across( where(is.numeric), .fns=sum )) %>% # sum up genewise feature counts
        ungroup()
      
      # From here we need to cast the motifcounts for the promoterregions, so that in the end we have one row per gene (instead of 1-2)
      motifdf_aggr <- motifdf_aggr %>% tidyr::pivot_wider(id_cols=c(label,seqnames,anno), 
                                                          names_from=class, 
                                                          values_from = !c(label,seqnames,anno, class),
                                                          values_fill = 0)
    } else {
      motifdf_aggr <-
        motifdf %>% 
        dplyr::select(!c(unselect_col,"name","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","gene_biotype","mgi_symbol")) %>%
        { if(weightby!=FALSE) mutate(., across(where(is.numeric), ~ (.x * get(weightby)))) else .} %>% # weight features by score
        dplyr::select(!any_of(as.character(weightby))) %>% # then we can drop the score since it will be nonsensical after the aggregation anyways
        group_by(label,seqnames,anno) %>%
        dplyr::summarise(across( where(is.numeric), .fns=sum )) %>% # sum up genewise feature counts
        ungroup()
    }
    
    
  } else { # loop for prox based assignments
    motifdf_aggr <-
      motifdf %>% 
      dplyr::select(!c("name","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","gene_biotype")) %>% 
      group_by(label,seqnames,anno) %>%
      dplyr::summarise(across( where(is.numeric), .fns=sum )) %>% # sum up genewise feature counts
      ungroup()
  }
  
  return(motifdf_aggr)
  
}


#------------function needs------------:
# featurematrix and indexes to perform split of test vs trainval set
# (trainval set is used for CV to find optimal regularization)
#------------function will------------:
# run cv.glmnet (as elastic net regression) and pick lambda.1se as regularization
# determine model performance on test set
# plots for raw counts of certain factors are currently turned of
complete_GLM_analysis <- function(motifdf, trainvalidx, genenames){
  
  # Create training subset for model development & testing set for model performance testing
  # make it based on chromosomes, instead of completely random
  # as in bpnet, we set aside chr 1,8 and 9 for testing
  
  #inTrain <- sort(sample(nrow(motifdata_abc_aggr), nrow(motifdata_abc_aggr)*0.75))
  
  features_train <- motifdf[ trainvalidx, -c(1,2,3)] %>% as.matrix()
  features_test <- motifdf[ -trainvalidx, -c(1,2,3)] %>% as.matrix()
  targets_train <- motifdf[ trainvalidx, ] %>% pull(label) %>% as.numeric(levels(.))[.] %>% as.matrix()
  targets_test <- motifdf[ -trainvalidx, ] %>% pull(label) %>% as.numeric(levels(.))[.] %>% as.matrix()
  
  #------------------------
  ## Elastic net
  #------------------------
  
  cvfit_net <- glmnet::cv.glmnet(x=features_train, 
                                 y=targets_train, 
                                 family="binomial", 
                                 type.measure = "auc",
                                 nfolds = 6,
                                 alpha=0.5)
  
  
  # Performance on training data
  targets_train_net.prob <- predict(cvfit_net,
                                    type="response",
                                    newx = features_train,
                                    s = 'lambda.min')
  pred_train_net <- ROCR::prediction(targets_train_net.prob[,1], targets_train) #only need the probabilities for 1's
  
  auc_ROCR_train_net <- ROCR::performance(pred_train_net, measure = "auc") #to assess AUC for model
  auc_ROCR_train_net@y.values[[1]]
  
  # Predict on unseen data
  targets_net.prob <- predict(cvfit_net,
                              type="response",
                              newx = features_test,
                              s = 'lambda.min')
  pred_net <- ROCR::prediction(targets_net.prob[,1], targets_test) #only need the probabilities for 1's
  
  auc_ROCR_net <- ROCR::performance(pred_net, measure = "auc") #to assess AUC for model
  auc_ROCR_net@y.values[[1]]
  
  #------------------------
  ## Model coefficients
  #------------------------
  
  coef_df_net <- coerce_coef2df( coefficients(cvfit_net, s = 'lambda.min'))
  
  gg_net <- ggplot( data=coef_df_net %>% filter(abs(estimates)>0) )+
    geom_bar(aes(x=reorder(names, -estimates), y=estimates),
             stat="identity")+
    labs(title=" ",x="", y="Net coefs")+
    #coord_flip()+
    theme(axis.text.x = element_text(angle = 45, hjust=1),
          plot.margin = unit(c(0, 0, -6, 10), "points"))
  
  #------------------------
  # ROC curves
  #------------------------
  # calculate probabilities for TPR/FPR for predictions
  perf_train_net <- ROCR::performance(pred_train_net,"tpr","fpr")
  perf_net <- ROCR::performance(pred_net,"tpr","fpr")
  
  # plot ROC curve
  
  gg_ROC <- ggplot()+
    geom_line( aes(x=perf_train_net@x.values[[1]], y=perf_train_net@y.values[[1]], colour = "train_net") ) +
    geom_line( aes(x=perf_net@x.values[[1]],y=perf_net@y.values[[1]], colour = "test_net") ) +
    geom_abline(intercept=0,slope=1, linetype=4, colour="grey")+
    scale_colour_manual(values=c("darkblue", "blue"), 
                        name=" ",
                        breaks=c("train_net", "test_net"),
                        labels=c(paste("net train. AUC:",round( auc_ROCR_train_net@y.values[[1]], 2)) ,
                                 paste("net test. AUC:", round( auc_ROCR_net@y.values[[1]], 2))
                        )
    )+
    labs(x="False positive rate", y="True positive rate")+
    theme(
      legend.position=c(0.65,0.15)
    )
  
  # add metric to global variable
  new_metrics = data.frame(
    net_train = auc_ROCR_train_net@y.values[[1]],
    net_test = auc_ROCR_net@y.values[[1]] 
  )
  
  
  full_panel <- ggpubr::ggarrange(gg_ROC,gg_net, nrow=2,
                                  labels=c("A","B"),
                                  heights = c(1,1))
  
  results <- list()
  results[[1]] <- full_panel
  results[[2]] <- new_metrics
  results[[3]] <- coef_df_net
  results[[4]] <- cvfit_net
  return(results)
}



#------------------------------------------------
## initialize object to track performance metrics
#------------------------------------------------

AUC_metrics = data.frame(
  # parameter combination
  condition=character(),
  motifdata=character(),
  excludepromoters=character(),
  onlymax=character(),
  weight=character(),
  sepPromEnh=character(),
  # performance
  net_train = numeric(),
  net_test = numeric()
)

net_model_coefs = list()


#------------------------------------------------
## run proximity based
#------------------------------------------------

# This is independent of the ABC results (no need to loop through different assignment variations)

print("Running GLM on prox-based assignments")
my_rdsfile <- paste0(opt$outdir,"prox.rds")
if(!file.exists(my_rdsfile)){
  
  motifdata_aggr <- merge_motifdata_with_assignments(motifcounts_summitregion,
                                                     assignment_summit_prox,
                                                     contrast_DexLPSvLPS,
                                                     maxonly="prox",
                                                     excludepromoters=FALSE,
                                                     weightby = FALSE,
                                                     sepPromEnh = FALSE)
  motifdata_aggr_scaled <- motifdata_aggr %>% mutate(., across(where(is.numeric), ~(scale(.) %>% as.vector)))
  
  motifdata_aggr_tranval_idx <- motifdata_aggr_scaled %>% with(which(seqnames!="chr1" & seqnames!="chr8" & seqnames!="chr9"))
  
  performance <- complete_GLM_analysis (motifdata_aggr_scaled, motifdata_aggr_tranval_idx, genenames=genenames)
  
  #------------------------
  ## Raw counts
  #------------------------
  nr3c1counts <- motifdata_aggr %>% 
    group_by(label) %>% 
    dplyr::count(NR3C1) %>% 
    mutate(counts_fac = case_when(NR3C1>=1 ~ ">=1",
                                  NR3C1==0 ~ "0"), # aggregate it towards to top end
           counts_fac = factor(as.character(counts_fac), levels=c("0",">=1")) ) %>%
    group_by(counts_fac) %>% 
    mutate(freq = n / sum(n)) %>%
    group_by(label,counts_fac) %>%
    dplyr::summarize(sum_n=sum(n),
              sum_freq=sum(freq))
  
  gg_nr3c1counts <- ggplot(data=nr3c1counts , aes(x=counts_fac,y=sum_freq,fill=label)) + 
    geom_bar(stat="identity", position="dodge", alpha=0.5)+
    geom_text(aes(label=format(sum_freq, digits=2)),
              position = position_dodge(.9), 
              vjust = -0.5, 
              size = 3) + 
    scale_fill_manual(values=c("0"="blue", "1"="red"),
                      labels=c("DOWN","UP"))+
    labs(x="NR3C1 matches", y="% genes", fill="")+
    theme(
      plot.margin = unit(c(10, 6, 0, 10), "points")
    )
  gg_nr3c1counts
  
  relcounts <- motifdata_aggr %>% 
    group_by(label) %>% 
    dplyr::count(REL) %>% 
    mutate(counts_fac = case_when(REL>=1 ~ ">=1",
                                  REL==0 ~ "0"), # aggregate it towards to top end
           counts_fac = factor(as.character(counts_fac), levels=c("0",">=1")) ) %>%
    group_by(counts_fac) %>%
    mutate(freq = n / sum(n)) %>%
    group_by(label,counts_fac) %>%
    dplyr::summarize(sum_n=sum(n),
              sum_freq=sum(freq))
  
  gg_relcounts <- ggplot(data=relcounts , aes(x=counts_fac,y=sum_freq,fill=label)) + 
    geom_bar(stat="identity", position="dodge", alpha=0.5)+
    geom_text(aes(label=format(sum_freq, digits=2)),
              position = position_dodge(.9), 
              vjust = -0.5, 
              size = 3) + 
    scale_fill_manual(values=c("0"="blue", "1"="red"),
                      labels=c("DOWN","UP"))+
    labs(x="REL matches", y="% genes", fill="")+
    theme(
      plot.margin = unit(c(10, 6, 0, 10), "points")
    )
  gg_relcounts
  
  gg_rawcounts <-  ggpubr::ggarrange( gg_nr3c1counts, gg_relcounts, ncol=2, common.legend = TRUE, legend="bottom")
  saveRDS(gg_rawcounts, paste0(opt$outdir,"gg_rawcounts.rds"))
  
  #------------------------

  
  saveRDS(performance,
          file=my_rdsfile)
} else {
  performance <- readRDS(my_rdsfile)
}

plot(performance[[1]])

AUC_metrics <- rbind(AUC_metrics,
                     "prox" = c(
                       condition="dexlps",
                       motifdata="motifcounts_summitregion",
                       excludepromoters=FALSE,
                       onlymax="prox",
                       weight=FALSE,
                       sepPromEnh=FALSE,
                       performance[[2]]
                       )
                     )

net_model_coefs[["prox"]] <- performance[[3]]


## Raw counts
gg_firstgene <- ggplot() + 
  geom_bar(data=motifdata_aggr , aes(NR3C1,fill=label), position="dodge", alpha=0.5)+
  scale_fill_manual(values=c("0"="blue", "1"="red"),
                    labels=c("DOWN","UP"))+
  labs(x="NR3C1", y="counts", fill="")+
  theme(
    plot.margin = unit(c(10, 6, 0, 10), "points")
  )
gg_firstgene

#------------------------------------------------
## run example
#------------------------------------------------


# motifdata_aggr <- merge_motifdata_with_assignments(
#   motifcounts_summitregion,
#   assignment_summitregion_abc_dexlps,
#   contrast_DexLPSvLPS,
#   maxonly=FALSE, 
#   excludepromoters=FALSE, 
#   weightby=FALSE, 
#   sepPromEnh=FALSE)
# 
# motifdata_aggr_tranval_idx <- motifdata_aggr %>% with(which(seqnames!="chr1" & seqnames!="chr8" & seqnames!="chr9"))
# performance <- complete_GLM_analysis (motifdata_aggr, motifdata_aggr_tranval_idx, genenames=genenames)
# plot(performance[[1]])
# AUC_metrics <- rbind(AUC_metrics,
#                      "test" = c(condition="dexlps",
#                                 motifdata="motifcounts_peak",
#                                 onlymax=FALSE,
#                                 excludepromoters=FALSE,
#                                 weight=FALSE,
#                                 sepPromEnh=FALSE,
#                                 performance[[2]]))
# net_model_coefs[["test"]] <- performance[[3]]
# 
# gg_firstgene <- ggplot() + 
#   geom_bar(data=performance[[5]] , aes(NR3C1,fill=label), position="dodge", alpha=0.5)+
#   scale_fill_manual(values=c("0"="blue", "1"="red"),
#                     labels=c("DOWN","UP"))+
#   labs(x="NR3C1", y="counts", fill="")+
#   theme(
#     plot.margin = unit(c(10, 6, 0, 10), "points")
#   )
# gg_firstgene


#------------------------------------------------
## peakregion_counts - proxbased assignment
## peakregion_counts - abcbased assignment (both conditions)
## enhancerregion_counts - abcbased assignment (2 conditions)
#------------------------------------------------

# rewrite to assign values for both conditions and then use those to compute the difference as well

not_all_na <- function(x) any(!is.na(x))

for (motifdata in c("motifcounts_abcregion","motifcounts_summitregion")){
  for(excludepromoters in c(FALSE,"all","onlyNONself")){
    for (onlymax in c(TRUE,FALSE)){
      for (sepPromEnh in c(TRUE,FALSE)){
        for (weight in c(FALSE, "abcscore")){
          
          # set assignments fitting for the input data
          if(motifdata=="motifcounts_abcregion"){
            assignment_dexlps <- assignment_abcregion_dexlps
            assignment_lps <- assignment_abcregion_lps
            motifcounts_dexlps <- motifcounts_abcregion_dexlps
            motifcounts_lps <- motifcounts_abcregion_lps
          } else if (motifdata=="motifcounts_summitregion" ){
            # in this case the motifcounts for the 2 conditions are the same, but there assignments differ
            assignment_dexlps <- assignment_summits_abcregion_dexlps
            assignment_lps <- assignment_summits_abcregion_lps
            motifcounts_dexlps <- motifcounts_summitregion
            motifcounts_lps <- motifcounts_summitregion
          } else {break}
          
          #-----------------------DEXLPS------------------------------
          featurematrix_dexlps <- merge_motifdata_with_assignments(motifcounts_dexlps,
                                                                   assignment_dexlps,
                                                                   contrast_DexLPSvLPS,
                                                                   weightby=weight, 
                                                                   maxonly=onlymax, 
                                                                   excludepromoters=excludepromoters, 
                                                                   sepPromEnh=sepPromEnh)
          
          # remove columns that are NA after scaling 
          featurematrix_dexlps_scaled <- featurematrix_dexlps %>% 
            mutate(., across(where(is.numeric), ~(scale(.) %>% as.vector))) %>%
            dplyr::select(where(not_all_na))
          
          #-----
          
          motifdata_aggr_tranval_idx <- featurematrix_dexlps_scaled %>% with(which(seqnames!="chr1" & seqnames!="chr8" & seqnames!="chr9"))
          modelname <- paste(motifdata,"condition_dexlps_exclprom",excludepromoters,"onlymax",onlymax,"sepPromEnh",sepPromEnh,"weight",weight, sep="_")
          print(modelname)
          
          my_rdsfile <- here( paste0(opt$outdir, modelname,".rds"))
          if(!file.exists(my_rdsfile)){
            performance <- complete_GLM_analysis (featurematrix_dexlps_scaled, motifdata_aggr_tranval_idx, genenames=genenames)
            saveRDS(performance,file=my_rdsfile)
          } else {
            performance <- readRDS(my_rdsfile)
          }
          
          AUC_metrics <- rbind(AUC_metrics,
                               c(
                                 condition="dexlps",
                                 motifdata=motifdata,
                                 excludepromoters=excludepromoters,
                                 onlymax=onlymax,
                                 weight=weight,
                                 sepPromEnh=sepPromEnh,
                                 performance[[2]]
                                 )
                               )%>%
            magrittr::set_rownames(c(rownames(AUC_metrics),modelname))
          
          net_model_coefs[[modelname]] <- performance[[3]]
          
          #--------------------------LPS---------------------------
          featurematrix_lps <- merge_motifdata_with_assignments(motifcounts_lps,
                                                                assignment_lps,
                                                                contrast_DexLPSvLPS,
                                                                weightby=weight, 
                                                                maxonly=onlymax, 
                                                                excludepromoters=excludepromoters, 
                                                                sepPromEnh=sepPromEnh)
          featurematrix_lps_scaled <- featurematrix_lps %>% 
            mutate(., across(where(is.numeric), ~(scale(.) %>% as.vector))) %>%
            dplyr::select(where(not_all_na))
          #-----
          motifdata_aggr_tranval_idx <- featurematrix_lps_scaled %>% with(which(seqnames!="chr1" & seqnames!="chr8" & seqnames!="chr9"))
          modelname <- paste(motifdata,"condition_lps_exclprom",excludepromoters,"onlymax",onlymax,"sepPromEnh",sepPromEnh,"weight",weight, sep="_")
          print(modelname)
          
          my_rdsfile <- here(paste0(opt$outdir, modelname,".rds"))
          if(!file.exists(my_rdsfile)){
            performance <- complete_GLM_analysis (featurematrix_lps_scaled, motifdata_aggr_tranval_idx, genenames=genenames)
            saveRDS(performance,file=my_rdsfile)
          } else {
            performance <- readRDS(my_rdsfile)
          }
          
          AUC_metrics <- rbind(AUC_metrics,
                               c(
                                 condition="lps",
                                 motifdata=motifdata,
                                 excludepromoters=excludepromoters,
                                 onlymax=onlymax,
                                 weight=weight,
                                 sepPromEnh=sepPromEnh,
                                 performance[[2]]
                                 )
                               )%>%
            magrittr::set_rownames(c(rownames(AUC_metrics),modelname))
          
          net_model_coefs[[modelname]] <- performance[[3]]
          
          
          #-----------------------DIFFERENCE------------------------------
          
          # they contain the same motifs, but not the same genes. 
          # doublecheck that all columnnames are identical
          table(colnames(featurematrix_dexlps) == colnames(featurematrix_lps))
          
          #  motifcounts missing in one condition should be set to 0 
          
          merged_featurematrix <- merge(featurematrix_dexlps,
                                        featurematrix_lps, 
                                        by=c("anno","label","seqnames"),
                                        all=TRUE)
          
          # replace missing counts in one of the conditions with 0
          merged_featurematrix[is.na(merged_featurematrix)] <- 0
          
          # use dataframe suffix to grab respective columns
          featurematrix_diff <- 
            merged_featurematrix[grep(".x$",colnames(merged_featurematrix))] - 
            merged_featurematrix[grep(".y$",colnames(merged_featurematrix))] 
          # tidy up column names
          colnames(featurematrix_diff) <- gsub(".x$","", colnames(featurematrix_diff))
          
          # add first 3 columns back after computing the difference of the counts
          featurematrix_diff <- cbind(merged_featurematrix[1:3],featurematrix_diff)
          
          featurematrix_diff_scaled <- featurematrix_diff %>% 
            mutate(., across(where(is.numeric), ~(scale(.) %>% as.vector))) %>%
            dplyr::select(where(not_all_na))
          
          # run GLM and look at performance  
          featurematrix_diff_tranval_idx <- featurematrix_diff_scaled %>% with(which(seqnames!="chr1" & seqnames!="chr8" & seqnames!="chr9"))
          
          modelname <- paste(motifdata,"condition_DexLPS-LPS_exclprom",excludepromoters,"onlymax",onlymax,"sepPromEnh",sepPromEnh,"weight",weight, sep="_")
          print(modelname)
          my_rdsfile <- here(paste0(opt$outdir, modelname,".rds"))
          if(!file.exists(my_rdsfile)){
            performance <- complete_GLM_analysis (featurematrix_diff_scaled, featurematrix_diff_tranval_idx, genenames=genenames)
            saveRDS(performance,file=my_rdsfile)
          } else {
            performance <- readRDS(my_rdsfile)
          }
          
          AUC_metrics <- rbind(AUC_metrics,
                               c(
                                 condition="DexLPS-LPS",
                                 motifdata=motifdata,
                                 excludepromoters=excludepromoters,
                                 onlymax=onlymax,
                                 weight=weight,
                                 sepPromEnh=sepPromEnh,
                                 performance[[2]]
                                 )
                               )%>%
            magrittr::set_rownames(c(rownames(AUC_metrics),modelname))
          
          net_model_coefs[[modelname]] <- performance[[3]]
          
          
        }
      }
    }
  }
}


#------------------------------------------------
## put model coefs into dataframes
#------------------------------------------------
motifnames <- colnames(motifcounts_summitregion)[3:ncol(motifcounts_summitregion)]
  
# make empty dataframe with all motifs
model_coefs_joint <- data.frame(featurename=motifnames)

# when we treat promoters and enhancers separately, the featurenames differ
model_coefs_sep <- data.frame(featurename= c(paste(motifnames,"genic",sep="_"),
                                             paste(motifnames,"promoter",sep="_"),
                                             paste(motifnames,"intergenic",sep="_"))
)

for (model in names(net_model_coefs)){
  if (grepl("sepPromEnh_FALSE",model)){
    model_coefs_joint <- merge(model_coefs_joint,
                               net_model_coefs[[model]],
                               by.x="featurename", by.y="names",
                               all=TRUE) %>% 
      dplyr::rename(!! model := "estimates")
  } else if (model=="prox"){
    model_coefs_joint <- merge(model_coefs_joint,
                               net_model_coefs[[model]],
                               by.x="featurename", by.y="names",
                               all=TRUE) %>% 
      dplyr::rename(!! model := "estimates")
  } else if (grepl("sepPromEnh_TRUE",model)){
    model_coefs_sep <- merge(model_coefs_sep,
                             net_model_coefs[[model]],
                             by.x="featurename", by.y="names",
                             all=TRUE) %>% 
      dplyr::rename(!! model := "estimates")
  } else { stop("The modelnames don't match the expected pattern.") }
}

#------------------------------------------------
## save data for plotting
#------------------------------------------------
saveRDS (model_coefs_joint, paste0(opt$outdir, "model_coefs_joint.rds"))
saveRDS (model_coefs_sep, paste0(opt$outdir, "model_coefs_sep.rds"))
saveRDS (AUC_metrics, paste0(opt$outdir, "AUC_metrics.rds"))

#------------------------------------------------
sessionInfo()
