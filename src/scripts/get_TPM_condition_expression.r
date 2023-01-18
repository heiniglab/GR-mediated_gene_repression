suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))

option_list <- list(
  make_option(c("-g", "--genekey"),
              type="character",
              help="Path to biomart derived genekey"),
  make_option(c("-e", "--expression"),
              type="character",
              help="Path to TPM expression matrix")
)

suppressPackageStartupMessages(library(dplyr, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts=F, quietly=T))

opt <- parse_args(OptionParser(option_list=option_list))

tpm_counts <- as.matrix(read.table(
  opt$expression,
  quote = "\"", 
  header=TRUE, 
  row.names = 1 ))


geneKey_GRCm38.p6 <- read.table(
  opt$genekey,
  header=TRUE,
  sep="\t")

tpm_counts_ext <- merge(as.data.frame(tpm_counts),
                        geneKey_GRCm38.p6, 
                        by.x="row.names", 
                        by.y="ensembl_gene_id", 
                        all.x=TRUE)


# compute the mean per condition
tpm_counts_ext <-
  tpm_counts_ext %>% 
  rowwise() %>%
  mutate(mean_dexlps=mean(c(LPS_Dex1_nascent, LPS_Dex2_nascent, LPS_Dex3_nascent))) %>%
  mutate(mean_lps=mean(c(LPS1_nascent, LPS2_nascent, LPS3_nascent ))) %>%
  mutate(mean_veh=mean(c(V1_nascent, V2_nascent )))


write.table(tpm_counts_ext[,c("Row.names","mean_veh")],
            file="results/current/abcmodel/expression/Veh_tpm.tsv",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            sep="\t")

write.table(tpm_counts_ext[,c("Row.names","mean_lps")],
            file="results/current/abcmodel/expression/LPS_tpm.tsv",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            sep="\t")

write.table(tpm_counts_ext[,c("Row.names","mean_dexlps")],
            file="results/current/abcmodel/expression/DexLPS_tpm.tsv",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            sep="\t")


