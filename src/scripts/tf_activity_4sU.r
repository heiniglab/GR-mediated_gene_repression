suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))

option_list <- list(
  make_option(c("--gene_annot"),
              type="character",
              help="Path to gencode annotation used to retrieve promoter regions"),
  make_option(c("--binding_sites_remap"),
              type="character",
              help="Path to peaks from the 2022 mouse remap release filtered for macrophages"),
  make_option(c("--rna_nascent_fpkm"),
              type="character",
              help="Path to normalized expression counts of nascent samples"),
  make_option(c("--genekey"),
              type="character",
              help="Path to genekey used to map ensembl geneIDs to mgi symbols"),
  make_option(c("--cage"),
              type="character",
              help="Path to bed file with location of max score within each cage read cluster"),
  make_option(c("--outdir"),
              type="character",
              help="Path to output directory")
)

opt <- parse_args(OptionParser(option_list=option_list))

dir.create( opt$outdir )
# created the first time executing
ftfbs_tss_annot_rds <- paste0( opt$outdir, "tfbs_tss2000u200d_annot.rds")
ftfbs_tssCAGE_annot_rds <- paste0( opt$outdir, "tfbs_tssCAGE2000u200d_annot.rds")

suppressPackageStartupMessages(library(here, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(dplyr, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(plsgenomics, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(data.table, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(GenomicRanges, warn.conflicts=F, quietly=T))

# ------------------------------------------------------------------------------
print("Get TFBS to TSS assignments.")
# ------------------------------------------------------------------------------

# MAC specific TSS
mac_spec_tss <- rtracklayer::import.bed(opt$cage)
mac_spec_tss <- mac_spec_tss %>% 
  filter(!duplicated(mac_spec_tss$SYMBOL)) %>%
  plyranges::as_granges()

load_gene_annotation <- function(gene_annot) {

  # load gene annotation
  ga <- fread(gene_annot)
  
  # file format is: chr origin type start stop U strand U add_info
  colnames(ga) <- c("chr", "origin", "type", "start", "stop", "score", "strand",
                    "frame", "info")
  # extract ranges
  ra <- with(ga, GRanges(chr, IRanges(start, stop), strand))
  
  # extract the additional attributes and merge with ranges object
  attrs <- strsplit(ga$info, ";")
  gene_id <- sapply(attrs, function(x) { sapply(strsplit(x[grepl("gene_id",x)], " "), "[[", 2) })
  gene_name <- sapply(attrs, function(x) { sapply(strsplit(x[grepl("gene_name",x)], " "), "[[", 3) })
  gene_biotype <- sapply(attrs, function(x) { sapply(strsplit(x[grepl("gene_type",x)], " "), "[[", 3) })
  
  # remove any lingering quotes
  gene_id <- gsub("\"", "", gene_id)
  gene_name <- gsub("\"", "", gene_name)
  gene_biotype <- gsub("\"", "", gene_biotype)
  
  # add to ranges object
  names(ra) <- gene_id
  ra$SYMBOL <- gene_name
  ra$BIOTYPE <- gene_biotype
  
  # finally, filter out 'misc_RNA' types and unusual chromosomes
  ra <- ra[ra$BIOTYPE != "misc_RNA"]
  ra <- keepStandardChromosomes(ra)
  
  return(ra)
}

annotate_tfbs_to_tss <- function(binding_sites_remap,
                                 tss) {
  
  # get the TFBS regions from remap
  tfbs = rtracklayer::import(binding_sites_remap)
  
  ann = t(matrix(unlist(strsplit(values(tfbs)[,"name"], ",", fixed=T)), nrow=3))
  
  colnames(ann) = c("geo_id", "TF", "condition")
  
  values(tfbs) = DataFrame(name=values(tfbs)[,"name"],
                           data.frame(ann, stringsAsFactors=F))
  
  
  # create an annotation matrix for the TSS
  chip = paste(values(tfbs)[,"TF"], values(tfbs)[,"condition"], sep=".")
  chip_exp = unique(chip)
  
  tfbs_ann = sapply(chip_exp, function(x) overlapsAny(tss,
                                                      tfbs[chip == x]))
  rownames(tfbs_ann) = names(tss)
  
  return(tfbs_ann)
}

# TSS only from ref
#------------------
ga <- load_gene_annotation(opt$gene_annot)
tss <- promoters(ga, 2000, 200)
names(tss) <- tss$SYMBOL
  
tfbs_annot <- annotate_tfbs_to_tss(opt$binding_sites_remap, tss)
saveRDS(tfbs_annot, file=ftfbs_tss_annot_rds)

# TSS from CAGE and only from ref where we have none
#------------------

ga_noCAGE <- ga %>% filter(!SYMBOL %in% mac_spec_tss$SYMBOL)
# add info on those genes where we don't have mac specific TSS
tss_wCAGE <- c(mac_spec_tss,ga_noCAGE)
tss_wCAGE <- promoters(tss_wCAGE, 2000, 200)
names(tss_wCAGE) <- tss_wCAGE$SYMBOL
tfbs_annot_CAGE <- annotate_tfbs_to_tss(opt$binding_sites_remap, tss_wCAGE)
saveRDS(tfbs_annot_CAGE, file=ftfbs_tssCAGE_annot_rds)

# ------------------------------------------------------------------------------
print("Read in expression data.")
# ------------------------------------------------------------------------------

expr_fpkm <- read.table(opt$rna_nascent_fpkm, header=TRUE)
expr_fpkm_nototal <- expr_fpkm [ !grepl('total|Total', colnames(expr_fpkm))]

# ------------------------------------------------------------------------------
print("Get ensembl to MGI mapping")
# ------------------------------------------------------------------------------

geneKey <- read.delim(opt$genekey)

# ------------------------------------------------------------------------------
print("Change gene annotation from Ensembl to MGI")
# ------------------------------------------------------------------------------

expr <- merge(expr_fpkm_nototal, geneKey[,c("ensembl_gene_id", "mgi_symbol")],
              by.x='Geneid', by.y='ensembl_gene_id')

expr <- expr %>% filter(mgi_symbol!="")

# summarize the expression values for those mgi symbols that have multiple entries
expr <- expr %>% 
  group_by(mgi_symbol) %>% 
  summarise(across(2:(ncol(expr)-1), mean)) %>%
  tibble::column_to_rownames("mgi_symbol")

write.table(expr, file=paste0( opt$outdir,"FPKMcounts_mgiaggr.tsv"),
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t",
            quote = FALSE)

# ------------------------------------------------------------------------------
print("Define annotation and data subsets.")
# ------------------------------------------------------------------------------
tfs <- unique(sapply(strsplit(colnames(tfbs_annot), "\\."), "[[", 1))

# one TF might have the same target measured more than once -> summarize
tss_annot_summarized <- sapply(tfs, function(tf) {
  rowSums(tfbs_annot[,grepl(paste0(tf, "\\."), colnames(tfbs_annot)), drop=F])
})

# get TFs and their targets (which targets are expressed)
targets <- intersect(rownames(tfbs_annot), rownames(expr))
# we skip the filtering step for TFS that are expressed, and work with all of them for now
# (partially because the genenames and TF protein names might not even match and we'll look at the results in more detail afterwards anyways)
#tf_sub <- tfs[tfs %in% rownames(expr)]

# get the annotation and data subsets
annot_sub <- tss_annot_summarized[targets,,drop=F]
#annot_sub <- annot_sub[, tf_sub]
data_sub <- expr[targets,]


# ------------------------------------------------------------------------------
print("Estimating TFAs using PLS/SIMPLS and substituting.")
# ------------------------------------------------------------------------------
TFA <- plsgenomics::TFA.estimate(annot_sub, data_sub)$TFA

rownames(TFA) <- colnames(annot_sub)
colnames(TFA) <- colnames(data_sub)


# ------------------------------------------------------------------------------
print("Saving results.")
# ------------------------------------------------------------------------------
write.table(file=paste0( opt$outdir, "TFA_4sU.tsv"), TFA,
            sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)

# ------------------------------------------------------------------------------
print("Same steps with CAGE specific data")
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
print("Define annotation and data subsets.")
# ------------------------------------------------------------------------------
tfs_CAGE <- unique(sapply(strsplit(colnames(tfbs_annot_CAGE), "\\."), "[[", 1))

# one TF might have the same target measured more than once -> summarize
tss_CAGE_annot_summarized <- sapply(tfs_CAGE, function(tf) {
  rowSums(tfbs_annot_CAGE[,grepl(paste0(tf, "\\."), colnames(tfbs_annot_CAGE)), drop=F])
})

# get TFs and their targets (which targets are expressed)
targets_CAGE <- intersect(rownames(tfbs_annot_CAGE), rownames(expr))
# we skip the filtering step for TFS that are expressed, and work with all of them for now
# (partially because the genenames and TF protein names might not even match and we'll look at the results in more detail afterwards anyways)
#tf_sub <- tfs[tfs %in% rownames(expr)]

# get the annotation and data subsets
annot_CAGE_sub <- tss_CAGE_annot_summarized[targets_CAGE,,drop=F]
#annot_sub <- annot_sub[, tf_sub]
data_sub_CAGE <- expr[targets_CAGE,]


# ------------------------------------------------------------------------------
print("Estimating TFAs using PLS/SIMPLS and substituting.")
# ------------------------------------------------------------------------------
TFA_CAGE <- plsgenomics::TFA.estimate(annot_CAGE_sub, data_sub_CAGE)$TFA

rownames(TFA_CAGE) <- colnames(annot_CAGE_sub)
colnames(TFA_CAGE) <- colnames(data_sub_CAGE)


# ------------------------------------------------------------------------------
print("Saving results.")
# ------------------------------------------------------------------------------
write.table(file=paste0( opt$outdir, "TFA_4sU_CAGE.tsv"), TFA_CAGE,
            sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)
