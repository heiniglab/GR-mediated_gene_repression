opt <- list()
opt$motifcounts_summitregion <- "results/current/GLM_dataprep/motifcounts_summitregion.rds"
opt$assignment_summit_prox <- "results/current/GLM_dataprep/assignment_summit_prox.rds"
opt$contrast_DexLPSvLPS <- "results/current/rnaseq_4sU/res_DexLPSvLPS_ext.tsv" 

suppressPackageStartupMessages(library(ggplot2, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(dplyr, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(here, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(ggcorrplot, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(RColorBrewer, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(circlize, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(ggpubr, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(ComplexHeatmap, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(clusterProfiler, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(DT, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(pathview, warn.conflicts=F, quietly=T))

set.seed(12345)

###########################
# Load data
###########################

motifcounts_summitregion <- readRDS ( opt$motifcounts_summitregion )
assignment_summit_prox <- readRDS( opt$assignment_summit_prox )
contrast_DexLPSvLPS <- read.delim(opt$contrast_DexLPSvLPS)


###########################
# Make featurematrix
###########################

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

motifdata_aggr_prox <- merge_motifdata_with_assignments(motifcounts_summitregion, assignment_summit_prox,contrast_DexLPSvLPS,maxonly="prox")

motifdata_aggr_prox_scaled <- motifdata_aggr_prox %>% mutate(., across(where(is.numeric), ~(scale(.) %>% as.vector)))

motifdata_aggr_prox_tranval_idx <- motifdata_aggr_prox_scaled %>% with(which(seqnames!="chr1" & seqnames!="chr8" & seqnames!="chr9"))

ncol( motifdata_aggr_prox_scaled[ , -c(1,2,3)]) # first three columns are labels, chromosome and gene annotation

targets_train <- motifdata_aggr_prox_scaled[ motifdata_aggr_prox_tranval_idx, ] %>% pull(label) %>% as.numeric(levels(.))[.] %>% as.matrix()
table(targets_train)[[0]]
targets_test <- motifdata_aggr_prox_scaled[ -motifdata_aggr_prox_tranval_idx, ] %>% pull(label) %>% as.numeric(levels(.))[.] %>% as.matrix()
table(targets_test)

###########################
# Pathway enrichment
###########################
#Check for pathway enrichment within geneset of genes that have any stathit

statcols <- colnames(motifdata_aggr_prox)[ grepl("Stat",colnames(motifdata_aggr_prox), ignore.case=TRUE) ]
statcols

motifdata_aggr_prox_stats <- motifdata_aggr_prox %>% select(c("anno",statcols)) %>% tibble::column_to_rownames(var="anno") %>% mutate(anystat_geneswithhit = rowSums(.)>0)


# filter for those that are downregulated
motifdata_aggr_prox_DOWN <- motifdata_aggr_prox%>% filter(label=="0")
# check correlation
cor(motifdata_aggr_prox_DOWN$`Stat5a::Stat5b`,motifdata_aggr_prox_DOWN$STAT3)
cor(motifdata_aggr_prox_DOWN$Stat5a,motifdata_aggr_prox_DOWN$STAT3)
cor(motifdata_aggr_prox_DOWN$Stat5b,motifdata_aggr_prox_DOWN$STAT3)

#--------------------

if (!file.exists(here("./results/current/kegg_enrichment/entrez_genekey.rds"))){
  my_ensembl <- biomaRt::useMart(biomart = 'ensembl',
                                 dataset="mmusculus_gene_ensembl",
                                 port=443)
  entrez_genekey <- biomaRt::getBM(mart=my_ensembl,
                                   attribute=c("ensembl_gene_id","gene_biotype","mgi_symbol", "entrezgene_id"))
  saveRDS(entrez_genekey,
          file=here("./results/current/kegg_enrichment/entrez_genekey.rds"))
} else {
  entrez_genekey <- readRDS(here("./results/current/kegg_enrichment/entrez_genekey.rds"))
}

mmu <- search_kegg_organism('mmu', by='kegg_code')

#entrez_geneids <- entrez_genekey[ entrez_genekey$ensembl_gene_id %in% (tt_interaction %>% filter(adj.P.Val<0.05) %>% dplyr::pull(geneid))  , "entrezgene_id" ]

motifdata_aggr_prox_stats <- merge(motifdata_aggr_prox_stats, entrez_genekey,
                           by.x="row.names",by.y="mgi_symbol")


kk_anystat <- enrichKEGG(gene = motifdata_aggr_prox_stats %>% filter(anystat_geneswithhit==TRUE) %>% dplyr::pull(entrezgene_id),
                             organism = 'mmu',
                             pvalueCutoff = 0.05)

as.data.frame(kk_anystat) %>% 
  dplyr::select(!geneID) %>% 
  mutate(across(where(is.numeric), ~format(.x, scientific=TRUE, digits=3))) %>% 
  relocate(p.adjust, .after="Description") %>% 
  head(n=50)


# visualize pathview with STAT3 counts
named_vector <- motifdata_aggr_prox_stats %>% filter(anystat_geneswithhit==TRUE) %>% dplyr::pull("STAT3")
names(named_vector) <- motifdata_aggr_prox_stats %>% filter(anystat_geneswithhit==TRUE) %>% dplyr::pull("entrezgene_id")


kegg_dir<-"/home/barbara.hoellbacher/projects/gr-manuscript/results/current/kegg_enrichment"
setwd(kegg_dir)
pathway_ID <- "mmu04061"
if(!file.exists(paste0("./results/current/kegg_enrichment/",pathway_ID,".pathview.png"))){
  pathview(gene.data  = named_vector,
           pathway.id = pathway_ID,
           species    = "mmu")
}

