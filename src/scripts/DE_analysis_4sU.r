
suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))

option_list <- list(
  make_option(c("-c", "--counts"),
              type="character",
              help="Path to count matrix"),
  make_option(c("-a", "--annotation"),
              type="character",
              help="Path to metadata with annotation"),
  make_option(c("-o", "--outdir"),
              type="character",
              help="Path to output directory"),
  make_option(c("-k", "--biomart_ensembl_version"),
              type="numeric",
              help="version of ensembl biomart to use for ensembl2mgi gene annotation"))
opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(library(dplyr, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(edgeR, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(DESeq2, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(biomaRt, warn.conflicts=F, quietly=T))

dir.create(opt$outdir)

#-------------------------------------------------
#---------Load data
#-------------------------------------------------
# Let's start by loading the count matrix holding the gene counts for all samples and by making a dataframe holding metadata such as treatment conditions.

counts <- as.matrix(read.table(opt$counts, quote = "\"", header=TRUE, row.names = 1 ))
anno <- read.delim(opt$annotation)

# parse relevant info from treatment variable
anno$treat <- anno$treatment %>%
  gsub(";4sU", "", .) %>%
  gsub("LPS, Dex", "LPS_Dex", .) %>%
  as.factor() %>%
  relevel(treat, ref = "Vehicle")

# parse total vs nascent from Foreign.ID
anno$group [ grepl('total|Total', anno$Foreign.ID) ] <- "total"
anno$group [ grepl('nascent', anno$Foreign.ID) ] <- "nascent"

# making sure anno is sorted to match column order in counts matrix
anno <- anno[ match(colnames(counts), 
                    anno$Foreign.ID),]
# check after resorting
all(anno$Foreign.ID == colnames(counts))

#-------------------------------------------------
#---------Filter genes
#-------------------------------------------------

# For now, we are not interested in making comparisons with the total RNA samples, so we exclude those before doing the normalization (and gene filtering) and running DESeq
anno <- anno[ anno$group == "nascent", ]

counts <- counts[,match(anno$Foreign.ID, colnames(counts))]

# do simple cpm computation without any library factors, to exclude lowly expressed genes
cpm_counts <- edgeR::cpm(counts, 
                         normalized.lib.sizes = FALSE)

# what does the read count distribution look like?
# compute median per gene and plot it as histogram
median_genecounts <- apply(cpm_counts, 1, FUN=median)

ggplot()+
  geom_histogram(aes(log(median_genecounts+1)), bins=30)

#Filter to keep genes that have a cpm of at least 0.2 in at least 1 samples

print("Filtering out lowly expressed genes")
keepRows <- rowSums(cpm_counts >=0.2) >= 1
table(keepRows)
median_genecounts <- median_genecounts[keepRows]

ggplot()+
  geom_histogram(aes(log(median_genecounts+1)), bins=30)

# use this to filter the gene counts matrix
counts <- counts[keepRows,]

#-------------------------------------------------
#---------Normalize counts
#-------------------------------------------------

# Create DESeq2Dataset object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = anno,
                              design= ~ treat)


# relevel to set Vehicle treated sample as reference
dds$treat <- relevel(dds$treat, ref = "Vehicle")

# this would also be done automatically when calling the DESeq() function
dds <- estimateSizeFactors(dds)

# NOTE: DESeq2 doesn’t actually use normalized counts, rather it uses the raw counts and models the normalization inside the Generalized Linear Model (GLM). 
# These normalized counts will be useful for downstream visualization of results, but cannot be used as input to DESeq2 or any other tools that perform differential expression analysis which use the negative binomial model.
normalized_counts <- counts(dds, normalized=TRUE)

write.table(normalized_counts, paste0(opt$outdir,"/DESeq_normalized_counts_nototal.tsv"), 
            col.names = TRUE,
            row.names = TRUE,
            sep="\t",
            quote=FALSE)

saveRDS(normalized_counts, paste0(opt$outdir,"/DESeq_normalized_counts_nototal.rds"))

write.table(anno, paste0(opt$outdir,"/anno_nototal.tsv"), 
            col.names = TRUE,
            row.names = FALSE,
            sep="\t",
            quote=FALSE)

#-------------------------------------------------
#---------Run DESeq2
#-------------------------------------------------
print("Running DESeq")
dds <- DESeq(dds)

res <- results(dds)

# based on the computed coefficients, we now compute a contrast of LPS vs LPS+Dex
res_DexLPSvLPS <- results(dds, contrast=c("treat","LPS_Dex","LPS"))
res_LPSvVeh <- results(dds, contrast=c("treat","LPS","Vehicle"))
res_DexLPSvVeh <- results(dds, contrast=c("treat","LPS_Dex","Vehicle"))

saveRDS(res_DexLPSvLPS, paste0(opt$outdir,file = "/contrast_DexVSDexLPS.rds"))
saveRDS(res_LPSvVeh, paste0(opt$outdir,file = "/contrast_LPSvVeh.rds"))
saveRDS(res_DexLPSvVeh, paste0(opt$outdir,file = "/contrast_DexLPSvVeh.rds"))

#-------------------------------------------------
#---------Annotate genes
#-------------------------------------------------

print("Annotating genes")

f_genekey <- paste0(opt$outdir,"/geneKey_biomart_mm_k",opt$biomart_ensembl_version,".txt")
if (!file.exists(f_genekey)){
  ensembl_mm <- biomaRt::useEnsembl(biomart = 'ensembl',
                                    dataset="mmusculus_gene_ensembl",
                                    version = as.numeric(opt$biomart_ensembl_version))
  # retrieve the geneKey to map ensembl IDs to mgi_symbols
  geneKey<- biomaRt::getBM(mart=ensembl_mm, attributes=c("ensembl_gene_id","gene_biotype","mgi_symbol"))
  write.table(geneKey, f_genekey, sep="\t")
  
} else {
  geneKey <- read.table(f_genekey, header=TRUE, sep="\t")
}
print("Done getting genekey")

#merge gene annotations to results table
res_DexLPSvLPS_ext <- merge(as.data.frame(res_DexLPSvLPS),
                            geneKey, 
                            by.x="row.names", 
                            by.y="ensembl_gene_id", 
                            all.x=TRUE)

res_LPSvVeh_ext <- merge(as.data.frame(res_LPSvVeh),
                         geneKey, 
                         by.x="row.names", 
                         by.y="ensembl_gene_id", 
                         all.x=TRUE)

res_DexLPSvVeh_ext <- merge(as.data.frame(res_DexLPSvVeh),
                            geneKey, 
                            by.x="row.names", 
                            by.y="ensembl_gene_id", 
                            all.x=TRUE)
#-------------------------------------------------
#---------Export contrasts
#-------------------------------------------------

print("Exporting the contrasts")

write.table(res_DexLPSvLPS_ext,paste0(opt$outdir,"/res_DexLPSvLPS_ext.tsv"),
            col.names = TRUE,
            row.names = FALSE,
            sep="\t",
            quote=FALSE)

write.table(res_LPSvVeh_ext,paste0(opt$outdir,"/res_LPSvVeh_ext.tsv"),
            col.names = TRUE,
            row.names = FALSE,
            sep="\t",
            quote=FALSE)

write.table(res_DexLPSvVeh_ext,paste0(opt$outdir,"/res_DexLPSvVeh_ext.tsv"),
            col.names = TRUE,
            row.names = FALSE,
            sep="\t",
            quote=FALSE)


