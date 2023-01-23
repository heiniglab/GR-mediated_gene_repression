suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))

option_list <- list(
  make_option(c("--auc"),
              type="character",
              help="Path to rds file with auc results of all models"),
  make_option(c("--dirname_featurematrizes"),
              type="character",
              help="directory that the unscaled featurematrizes were saved in"),
  make_option(c("--dirname_models"),
              type="character",
              help="directory that the trained models were saved in"),
  make_option(c("--outfig"),
              type="character",
              help="Path to output Figure")
)

opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(library(here, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(glmnet, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(dplyr, warn.conflicts=F, quietly=T))

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

# compare ROC curves for statistical significance
# which ones?
# let's start with the best performing one and our reference model

# what is the name of the best performing model?
AUC_metrics <- readRDS ( opt$auc ) 
bestmodel <- AUC_metrics %>% filter(net_test == max(AUC_metrics$net_test)) %>% rownames()
AUC_metrics[bestmodel,"net_train"]
# what is the performance of the best model on the training set?
#-----------------------------------------

# 0 write function that takes model and featurematrix as input and returns ROC object as output
not_all_na <- function(x) any(!is.na(x))
get_roc_object <- function(path_performance_rds, path_featurematrix_unscaled_rds){
  
  # 1. load cv.fit of model
  performance<- readRDS( path_performance_rds )
  cvfit_net <- performance[[4]]
  
  # 2. load testset of model
  motifdata <- readRDS( path_featurematrix_unscaled_rds )
  motifdata_scaled <- motifdata %>% 
    mutate(., across(where(is.numeric), ~(scale(.) %>% as.vector)))%>%
    dplyr::select(where(not_all_na))
  motifdata_tranval_idx <- motifdata_scaled %>% with(which(seqnames!="chr1" & seqnames!="chr8" & seqnames!="chr9"))
  
  features_train <- motifdata_scaled[ motifdata_tranval_idx, -c(1,2,3)] %>% as.matrix()
  features_test <- motifdata_scaled[ -motifdata_tranval_idx, -c(1,2,3)] %>% as.matrix()
  targets_train <- motifdata_scaled[ motifdata_tranval_idx, ] %>% pull(label) %>% as.numeric(levels(.))[.] %>% as.matrix()
  targets_test <- motifdata_scaled[ -motifdata_tranval_idx, ] %>% pull(label) %>% as.numeric(levels(.))[.] %>% as.matrix()
  
  # 3. Predict and return ROC
  targets_net.prob <- predict(cvfit_net,
                              type="response",
                              newx = features_test,
                              s = 'lambda.min')
  
  roc_object <- pROC::roc(as.factor(targets_test),targets_net.prob[,1], direction = "<")
  return(roc_object)
  
}

roc_reference <- get_roc_object ( paste0( opt$dirname_models, "prox.rds"),
                                  paste0( opt$dirname_featurematrizes, "prox.rds") )

roc_bestmodel <- get_roc_object ( paste0( opt$dirname_models, bestmodel,".rds"),
                                  paste0( opt$dirname_featurematrizes, bestmodel,".rds") )

# 4. Compare for pROC
# Should we compute the p value comparing our reference to all others and see how many others it beats?
#res <- pROC::roc.test( roc_reference, roc_bestmodel, method="delong", alternative="greater")
#res

#-------------------------------------------------------------------------------
#  comparison with reference model (proximity based)
#-------------------------------------------------------------------------------

# initialize a dataframe that will save the name of the model, the pvalue and whether the referene was better
all_pairwise_wprox <- data.frame(
  pvalues=numeric())

for (modelname in rownames(AUC_metrics) ) {
  # skip the prox model since this is what we compare to
  if (modelname=="prox"){next}
  print(modelname)
  roc_model <- get_roc_object ( paste0( opt$dirname_models, modelname,".rds"),
                                paste0( opt$dirname_featurematrizes, modelname,".rds") )
  roc_test_res <- pROC::roc.test( roc_reference, roc_model, method="delong", alternative="greater") 
    
  new_metrics <- data.frame(
    p.values=roc_test_res$p.value)
    
  all_pairwise_wprox <- rbind(all_pairwise_wprox,
                              new_metrics ) %>%
    magrittr::set_rownames(c(rownames(all_pairwise_wprox), modelname))
}

all_pairwise_wprox$input_is_summitregion <- grepl("motifcounts_summitregion*",rownames(all_pairwise_wprox))

# plot all p values to estimate where the histogram "levels out" and set this as lambda in Storey's q value
gg_storey_wprox <-
  ggplot(all_pairwise_wprox)+
  geom_histogram(aes(p.values,group=input_is_summitregion,fill=input_is_summitregion), bins=20 )+
  scale_fill_manual("",
                    labels = c("active regions", "GR summitregions"),
                    breaks = c("FALSE", "TRUE"),
                    values = c("chocolate4", "purple")) +
  scale_colour_manual("",
                      labels = c("active regions", "GR summitregions"),
                      breaks = c("FALSE", "TRUE"),
                      values = c("chocolate4", "purple")) +
  labs(x="p value",
       y="count")+
  geom_vline(xintercept=0.34) +
  theme(legend.position = c(0.7, 0.8))
gg_storey_wprox

storey_wprox <- qvalue::qvalue(all_pairwise_wprox$p.values,  lambda=0.34)
storey_wprox
print("Pi0 computed by Storey's method for comparison with reference model:")
print(storey_wprox$pi0)

table(all_pairwise_wprox$p.values<=0.05, all_pairwise_wprox$input_is_summitregion)

#-------------------------------------------------------------------------------
# comparison with best performing model
#-------------------------------------------------------------------------------

# initialize a dataframe that will save the name of the model, the pvalue and whether the referene was better
all_pairwise_wbest <- data.frame(pvalues=numeric())

for (modelname in rownames(AUC_metrics) ) {
  # skip the bestmodel since this is what we compare to
  if (modelname==bestmodel){
    print("Skipping self") 
    next}
  print(modelname)
  roc_model <- get_roc_object ( paste0(opt$dirname_models, modelname,".rds"),
                                paste0( opt$dirname_featurematrizes, modelname,".rds") )
  roc_test_res <- pROC::roc.test( roc_bestmodel, roc_model, method="delong", alternative="greater") 
  
  new_metrics <- data.frame(p.values=roc_test_res$p.value)
  
  all_pairwise_wbest <- rbind(all_pairwise_wbest,
                              new_metrics ) %>%
    magrittr::set_rownames(c(rownames(all_pairwise_wbest), modelname))
}

all_pairwise_wbest$input_is_summitregion <- grepl("motifcounts_summitregion*",rownames(all_pairwise_wbest))

# plot all p values to estimate where the histogram "levels out" and set this as lambda in Storey's q value
gg_storey_wbest <-
  ggplot(all_pairwise_wbest)+
  geom_histogram(aes(p.values,group=input_is_summitregion,fill=input_is_summitregion), bins=20 )+
  scale_fill_manual("",
                    labels = c("active regions", "GR summitregions"),
                    breaks = c("FALSE", "TRUE"),
                    values = c("chocolate4", "purple")) +
  scale_colour_manual("",
                      labels = c("active regions", "GR summitregions"),
                      breaks = c("FALSE", "TRUE"),
                      values = c("chocolate4", "purple")) +
  labs(x="p value",
       y="count")+
  geom_vline(xintercept=0.1) +
  theme(legend.position = c(0.7, 0.8))
gg_storey_wbest

storey_wbest <- qvalue::qvalue(all_pairwise_wbest$p.values,  lambda=0.1)
storey_wbest
print("Pi0 computed by Storey's method for comparison with best model:")
print(storey_wbest$pi0)

table(all_pairwise_wbest$p.values<=0.05)


#-------------------------------------------------------------------------------
# bivariate models with motifs of interest
#-------------------------------------------------------------------------------
# Load in featurematrix, scale it and fit model with out motifs of interest in bivariate model

motifdata <- readRDS( paste0(opt$dirname_featurematrizes, "prox.rds" ))
motifdata_scaled <- motifdata %>% 
  mutate(., across(where(is.numeric), ~(scale(.) %>% as.vector)))
motifdata_tranval_idx <- motifdata_scaled %>% with(which(seqnames!="chr1" & seqnames!="chr8" & seqnames!="chr9"))

features_train <- motifdata_scaled[ motifdata_tranval_idx, -c(1,2,3)] %>% as.matrix()
features_test <- motifdata_scaled[ -motifdata_tranval_idx, -c(1,2,3)] %>% as.matrix()
targets_train <- motifdata_scaled[ motifdata_tranval_idx, ] %>% pull(label) %>% as.numeric(levels(.))[.] %>% as.matrix()
targets_test <- motifdata_scaled[ -motifdata_tranval_idx, ] %>% pull(label) %>% as.numeric(levels(.))[.] %>% as.matrix()

model_coefficients <- data.frame(motifname = character(),
                                 intercept = numeric(),
                                 coefficient = numeric(),
                                 pvalue = numeric() )

for (motif in c("MEIS1","NFKB1","REL","POU2F1","TCF7","STAT3")){
  glm_res <- glm(targets_train ~ features_train[,c(motif)], family="binomial")
  new_model <- data.frame(motifname=motif,
                          intercept= glm_res$coefficients[1],
                          coefficient = glm_res$coefficients[2],
                          pvalue = coef(summary(glm_res))[2,4],
                          row.names = NULL)
  model_coefficients <- rbind( model_coefficients, new_model )
}


gg_bivariate <- ggplot(data=model_coefficients)+
  geom_point(aes(coefficient, -log10(pvalue)))+
  ggrepel::geom_label_repel(aes(coefficient, -log10(pvalue), 
                                label=motifname), size=2)+
  coord_cartesian(ylim=c(0,max(-log10(model_coefficients$pvalue))),
                  xlim=c(min(model_coefficients$coefficient),-0.1))+
  geom_hline(aes(yintercept=-log10(0.05)), linetype="dashed" )

#-------------------------------------------------------------------------------
# model performance on training set
#-------------------------------------------------------------------------------
melted_AUC_metrics <- AUC_metrics %>%
  reshape::melt(., measure.vars=c("net_train","net_test"))

#remove the proximity based model and add it in shape of a reference line instead
prox_reference_train <- melted_AUC_metrics %>% filter(variable=="net_train") %>% filter(onlymax=="prox") %>% pull(value) %>% as.numeric()
gg_AUC_train <- ggplot(melted_AUC_metrics %>% filter(variable=="net_train") %>% filter(onlymax!="prox"))+
  geom_bar(aes(x=weight, y=value, fill=condition, alpha=onlymax), 
           colour="black", size=0.5, width=0.9, 
           stat="identity",
           position=position_dodge())+
  geom_hline(aes(yintercept=prox_reference_train, linetype="reference"), colour="purple")+
  scale_linetype_manual(name="", values = c(2))+
  facet_grid( cols = vars(sepPromEnh,excludepromoters),
              rows= vars(motifdata),
              labeller = labeller(
                motifdata=c("motifcounts_abcregion"="active regions",
                            "motifcounts_summitregion"="GR summitregions"),
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

#-------------------------------------------------------------------------------
# arrange and save figure
#-------------------------------------------------------------------------------

gg_placeholder <- ggplot() + 
  theme(plot.margin = margin(t=0, l=0, r=0, b=0, unit = "cm"),
        axis.line = element_blank())

#A in this panel will be tha training performance

gg_c1 <- ggpubr::ggarrange(gg_storey_wbest, gg_storey_wprox, 
                           labels = c("B","C"),
                           ncol = 1, nrow = 2)
gg_c2 <- ggpubr::ggarrange(gg_bivariate, gg_placeholder, 
                           labels = c("D",NA),
                           ncol = 1, nrow = 2,
                           heights = c(1.5,1))
gg_r2 <-  ggpubr::ggarrange(gg_c1, gg_c2, 
                         labels = c(NA, NA),
                         ncol = 2, nrow = 1)

gg <-  ggpubr::ggarrange(gg_AUC_train, gg_r2, 
                         labels = c("A", NA),
                         ncol = 1, nrow = 2,heights = c(1,1))
gg
ggsave( opt$outfig, gg,
       width=190, height=200, units="mm",
       bg="white")
# also save it as pdf
ggsave( gsub(".png",".pdf",opt$outfig), gg,
        width=190, height=200, units="mm",
        bg="white")
