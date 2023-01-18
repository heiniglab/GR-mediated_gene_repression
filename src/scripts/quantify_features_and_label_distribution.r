# We didn't add this code to the "abc_prediction.r" script (even though code is largely shared) to maintain reproducibility of the code

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
  make_option(c("--model_coefs_joint"),
              type="character",
              help="Path to rds file with model coefficients of joint models"),
  make_option(c("--model_coefs_sep"),
              type="character",
              help="Path to rds file with model coefficients of models tha include enhancers and promoters separately"),
  make_option(c( "--featuredir"),
              type="character",
              help="Path to directory with unscales featurematrizes"),
  make_option(c( "--outfile"),
              type="character",
              help="Path to output file with metrics")
)

opt <- parse_args(OptionParser(option_list=option_list))

dir.create(opt$featuredir)

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

for (optname in names(opt)[2:11]){ #except for the featuredir
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

get_feature_and_label_metrics <- function(motifdf, trainvalidx, genenames){

  targets_train <- motifdf[ trainvalidx, ] %>% pull(label) %>% as.numeric(levels(.))[.] %>% as.matrix()
  targets_test <- motifdf[ -trainvalidx, ] %>% pull(label) %>% as.numeric(levels(.))[.] %>% as.matrix()
  
  metrics = data.frame(
    n_input_features = ncol( motifdf[ , -c(1,2,3)]),# first three columns are labels, chromosome and gene annotation
    n_neg_inst_train = table(targets_train)[['0']],
    n_pos_inst_train = table(targets_train)[['1']],
    n_neg_inst_test = table(targets_test)[['0']],
    n_pos_inst_test = table(targets_test)[['1']]
  )
  
  return(metrics)
}

#------------------------------------------------
## initialize object to gather all metrics
#------------------------------------------------
all_metrics=data.frame(
  n_input_features=numeric(),
  n_neg_inst_train=numeric(),
  n_pos_inst_train=numeric(),
  n_neg_inst_test=numeric(),
  n_pos_inst_test=numeric()
)
#------------------------------------------------
## run proximity based
#------------------------------------------------

#---------------------- 
# This is independent of the ABC results (no need to loop through different assignment variations)

print("Getting metrics of prox-based model")
my_rdsfile <- paste0(opt$featuredir,"prox.rds")
if(!file.exists(my_rdsfile)){
  
  motifdata_aggr_prox <- merge_motifdata_with_assignments(motifcounts_summitregion,
                                                     assignment_summit_prox,
                                                     contrast_DexLPSvLPS,
                                                     maxonly="prox",
                                                     excludepromoters=FALSE,
                                                     weightby = FALSE,
                                                     sepPromEnh = FALSE)
  
  saveRDS(motifdata_aggr_prox,file=my_rdsfile)
} else{
  motifdata_aggr_prox <- readRDS(my_rdsfile)
}

motifdata_aggr_prox_scaled <- motifdata_aggr_prox %>% mutate(., across(where(is.numeric), ~(scale(.) %>% as.vector)))
motifdata_aggr_prox_tranval_idx <- motifdata_aggr_prox_scaled %>% with(which(seqnames!="chr1" & seqnames!="chr8" & seqnames!="chr9"))
metrics_prox <- get_feature_and_label_metrics (motifdata_aggr_prox_scaled, motifdata_aggr_prox_tranval_idx, genenames=genenames)
 
all_metrics <- rbind(all_metrics,
                     "prox"=metrics_prox)


#---------------------- 

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
          modelname <- paste(motifdata,"condition_dexlps_exclprom",excludepromoters,"onlymax",onlymax,"sepPromEnh",sepPromEnh,"weight",weight, sep="_")
          print(modelname)
          
          my_rdsfile <- here( paste0(opt$featuredir, modelname,".rds"))
          if(!file.exists(my_rdsfile)){
            featurematrix_dexlps <- merge_motifdata_with_assignments(motifcounts_dexlps,
                                                                   assignment_dexlps,
                                                                   contrast_DexLPSvLPS,
                                                                   weightby=weight, 
                                                                   maxonly=onlymax, 
                                                                   excludepromoters=excludepromoters, 
                                                                   sepPromEnh=sepPromEnh)
            
            saveRDS(featurematrix_dexlps, my_rdsfile)
          } else{
            featurematrix_dexlps <- readRDS(my_rdsfile)
          }
          # remove columns that are NA after scaling
          featurematrix_dexlps_scaled <- featurematrix_dexlps %>%
            mutate(., across(where(is.numeric), ~(scale(.) %>% as.vector))) %>%
            dplyr::select(where(not_all_na))
          motifdata_aggr_tranval_idx <- featurematrix_dexlps_scaled %>% with(which(seqnames!="chr1" & seqnames!="chr8" & seqnames!="chr9"))
          new_metric <- get_feature_and_label_metrics (featurematrix_dexlps_scaled, motifdata_aggr_tranval_idx, genenames=genenames)
            
          all_metrics <- rbind(all_metrics,
                               new_metric)%>%
            magrittr::set_rownames(c(rownames(all_metrics),modelname))
          
          #--------------------------LPS---------------------------
          modelname <- paste(motifdata,"condition_lps_exclprom",excludepromoters,"onlymax",onlymax,"sepPromEnh",sepPromEnh,"weight",weight, sep="_")
          print(modelname)
          
          my_rdsfile <- here(paste0(opt$featuredir, modelname,".rds"))
          if(!file.exists(my_rdsfile)){
            featurematrix_lps <- merge_motifdata_with_assignments(motifcounts_lps,
                                                                assignment_lps,
                                                                contrast_DexLPSvLPS,
                                                                weightby=weight, 
                                                                maxonly=onlymax, 
                                                                excludepromoters=excludepromoters, 
                                                                sepPromEnh=sepPromEnh)

            saveRDS(featurematrix_lps,file=my_rdsfile)
          } else{
            featurematrix_lps <- readRDS(my_rdsfile)
          }
          featurematrix_lps_scaled <- featurematrix_lps %>% 
            mutate(., across(where(is.numeric), ~(scale(.) %>% as.vector))) %>%
            dplyr::select(where(not_all_na))
          motifdata_aggr_tranval_idx <- featurematrix_lps_scaled %>% with(which(seqnames!="chr1" & seqnames!="chr8" & seqnames!="chr9"))
          new_metric <- get_feature_and_label_metrics (featurematrix_lps_scaled, motifdata_aggr_tranval_idx, genenames=genenames)

          all_metrics <- rbind(all_metrics,
                               new_metric)%>%
            magrittr::set_rownames(c(rownames(all_metrics),modelname))
          
          #-----------------------DIFFERENCE------------------------------
          
          modelname <- paste(motifdata,"condition_DexLPS-LPS_exclprom",excludepromoters,"onlymax",onlymax,"sepPromEnh",sepPromEnh,"weight",weight, sep="_")
          print(modelname)
          my_rdsfile <- here(paste0(opt$featuredir, modelname,".rds"))
          if(!file.exists(my_rdsfile)){
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
            
            saveRDS(featurematrix_diff,file=my_rdsfile)
          } else{
            featurematrix_diff <- readRDS(my_rdsfile)
          }
          featurematrix_diff_scaled <- featurematrix_diff %>% 
            mutate(., across(where(is.numeric), ~(scale(.) %>% as.vector))) %>%
            dplyr::select(where(not_all_na))
          # run GLM and look at performance  
          featurematrix_diff_tranval_idx <- featurematrix_diff_scaled %>% with(which(seqnames!="chr1" & seqnames!="chr8" & seqnames!="chr9"))
          new_metric <- get_feature_and_label_metrics (featurematrix_diff_scaled, featurematrix_diff_tranval_idx, genenames=genenames)
          
          all_metrics <- rbind(all_metrics,
                               new_metric)%>%
            magrittr::set_rownames(c(rownames(all_metrics),modelname))
          
          
        }
      }
    }
  }
}

all_metrics <- all_metrics %>%
  mutate(ratio_pos_inst_train = n_pos_inst_train/(n_pos_inst_train+n_neg_inst_train),
         ratio_pos_inst_test = n_pos_inst_test/(n_pos_inst_test+n_neg_inst_test)
         )

# let'S add the number of non-zero coefficients to this
model_coefs_joint <- readRDS( opt$model_coefs_joint )
model_coefs_sep <- readRDS( opt$model_coefs_sep )

nonzero_columentries <- rbind ( colSums(model_coefs_joint != 0, na.rm=TRUE) %>% as.data.frame(),
                            colSums(model_coefs_sep != 0, na.rm=TRUE) %>% as.data.frame() )
colnames(nonzero_columentries) <- "n_sel_features_inclintercept"

all_metrics <- merge(all_metrics, nonzero_columentries, by="row.names") %>% 
  relocate(n_sel_features_inclintercept, .after = n_input_features )


# parse feature engineering choices from modelname and add them as variables
all_metrics <- all_metrics %>%
  mutate(condition = case_when(grepl("condition_lps", Row.names) ~ "LPS",
                              grepl("condition_DexLPS-LPS", Row.names) ~ "Dex+LPS - LPS",
                              grepl("condition_dexlps", Row.names) ~ "Dex+LPS"),
         input_region = case_when(grepl("motifcounts_summitregion", Row.names) ~ "GR summitregions",
                                  grepl("motifcounts_abcregion", Row.names) ~ "active regions"),
         excludepromoters = case_when(grepl("exclprom_all", Row.names) ~ "all",
                                      grepl("exclprom_onlyNONself", Row.names) ~ "nonself",
                                      grepl("exclprom_FALSE", Row.names) ~ "none"),
         mapping = case_when(grepl("onlymax_FALSE", Row.names) ~ "1-to-many (ABC-based)",
                             grepl("onlymax_TRUE", Row.names) ~ "1-to-1 (ABC-based)"),
         weight = case_when(grepl("weight_abcscore", Row.names) ~ "abcscore",
                             grepl("weight_FALSE", Row.names) ~ "none"),
         prom_and_enh_features = case_when(grepl("sepPromEnh_TRUE", Row.names) ~ "separate",
                                      grepl("sepPromEnh_FALSE", Row.names) ~ "aggregate")
         
         )
# manually update choices for the reference model
all_metrics <- all_metrics %>%
  rows_update( tibble(Row.names = "prox", 
                      condition="Dex+LPS",
                      input_region="GR summitregions",
                      mapping="1-to-1 (proximity-based)",
                      weight="none"), 
               by = "Row.names")

summary(all_metrics)

summary(all_metrics$ratio_pos_inst_train-all_metrics$ratio_pos_inst_test)

write.table(all_metrics,
            file=opt$outfile,
            quote=FALSE, sep="\t", row.names = FALSE, col.names = TRUE )
