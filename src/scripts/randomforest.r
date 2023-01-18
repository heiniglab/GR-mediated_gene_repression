library(randomForest)
library(dplyr)

set.seed(42)

opt <- list()
opt$outdir <- "results/current/GLM_featurematrizes_unscaled/"

# list with all "untuned" rf models
rf=list()

all_files <- list.files(opt$outdir, full.names = TRUE)
featurematrizes_path <- all_files[ - grepl("metrics_n_features_label_distro.tsv",all_files) ]

not_all_na <- function(x) any(!is.na(x))
#------------------------------------------------
## Loop through all featurematrizes
#------------------------------------------------

for (featurematrix_path in featurematrizes_path) {
  new_name <- basename(featurematrix_path)
  print(new_name)
  motifdata_aggr <- readRDS(featurematrix_path)
  motifdata_aggr_scaled <- motifdata_aggr %>% 
    mutate(., across(where(is.numeric), ~(scale(.) %>% as.vector))) %>%
    dplyr::select(where(not_all_na))
  motifdata_aggr_tranval_idx <- motifdata_aggr_scaled %>% with(which(seqnames!="chr1" & seqnames!="chr8" & seqnames!="chr9"))
  
  features_train <- motifdata_aggr_scaled[ motifdata_aggr_tranval_idx, -c(1,2,3)] %>% as.matrix()
  features_test <- motifdata_aggr_scaled[ -motifdata_aggr_tranval_idx, -c(1,2,3)] %>% as.matrix()
  targets_train <- motifdata_aggr_scaled[ motifdata_aggr_tranval_idx, ] %>% pull(label) %>% as.numeric(levels(.))[.] %>% as.matrix()
  targets_test <- motifdata_aggr_scaled[ -motifdata_aggr_tranval_idx, ] %>% pull(label) %>% as.numeric(levels(.))[.] %>% as.matrix()
  
  new_rf <- randomForest(y = factor(targets_train), 
                         x = features_train ,
                         ytest = factor(targets_test),
                         xtest = features_test
  )
  
  rf[[new_name]] <- new_rf
}

dir.create("results/current/randomforest")
saveRDS(rf, file="results/current/randomforest/rf.rds")

#------------------------------------------------
## Extract importance and performance matrics from models
#------------------------------------------------

#importance <- data.frame()

metrics <- data.frame(
  train_accuracy = numeric(),
  test_accuracy = numeric()
)

for(i in seq_along(1:length(names(rf)))) {
  name <- names(rf)[i]
  print(name)
  # merging doesn't work since some are joint and some are separate for prom and enh
  # MeanDecreaseGini
  #new_importance <- rf[[name]]$importance
  #colnames(new_importance) <- name
  #if(i==1){
  #  importance <- new_importance
  #} else{
  #  importance <- merge(importance,new_importance)
  #}
  
  cm_train <- rf[[name]]$confusion
  cm_test <- rf[[name]]$test$confusion
  
  new_metrics <- data.frame(
    train_accuracy = sum(cm_train[1], cm_train[4]) / sum(cm_train[1:4]),
    test_accuracy = sum(cm_test[1], cm_test[4]) / sum(cm_test[1:4])
  )
    
  rownames(new_metrics) <- name
  metrics <- rbind( metrics, new_metrics )
}
metrics
#View(importance)



#----------------------
# Hyperparameter tuning
# not finished
#--------------------

library(caret)

# Define the control
trControl <- trainControl(method = "cv",
                          number = 10,
                          search = "grid")

set.seed(1234)
# Run the model
rf_default <- train(y = factor(targets_train),
                    x = features_train,
                    method = "rf",
                    metric = "Accuracy",
                    trControl = trControl)
# Print the results
print(rf_default)

rf_default$finalModel


set.seed(1234)
tuneGrid <- expand.grid(.mtry = c(1: 50))
rf_mtry <- train(y = factor(targets_train),
                 x = features_train,
                 method = "rf",
                 metric = "Accuracy",
                 tuneGrid = tuneGrid,
                 trControl = trControl,
                 importance = TRUE)
print(rf_mtry)

