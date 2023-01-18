
suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))

option_list <- list(
  make_option(c("--log2fcthresh"),
              type="numeric",
              help="Log2FC threshold used in addition to adj.pval to define significant genes"),
  make_option(c("--chipseq_summits"),
              type="character",
              help="Path to summit file of IDR peaks"),
  make_option(c("--genekey"),
              type="character",
              help="Path to biomart genekey that mapps ensembl geeneIDs to MGI symbols"),
  make_option(c("--contrast_DexVSDexLPS"),
              type="character",
              help="Path to annotated tsv file of DeSeq2 contrast of DexLPS vs LPS"),
  make_option(c("--meme_db_path"),
              type="character",
              help="Path to JASPAR motif db file"),
  make_option(c( "--rna_nascent_fpkm"),
              type="character",
              help="FPKM matrix of 4sU experiment"),
  make_option(c("-o", "--outdir"),
              type="character",
              help="Path to output directory"))

opt <- parse_args(OptionParser(option_list=option_list))

# set output for logfile to retrieve stats for plot later
sink(file=paste0(opt$outdir,"figure_proxanno_prepdata.out"))

suppressPackageStartupMessages(library(memes, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(universalmotif, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(biomaRt, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(TxDb.Mmusculus.UCSC.mm10.knownGene, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(BSgenome.Mmusculus.UCSC.mm10, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(dplyr, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(ChIPseeker, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(stringr, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(plyranges, warn.conflicts=F, quietly=T))

#-------------------------------
## Import references
#-------------------------------
# for gene annotation
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
# for the sequence
# either use masked or unmasked (the mask does NOT seem to be for repeats though)
mm.genome <- BSgenome.Mmusculus.UCSC.mm10

#-------------------------------
### Determine expressed genes
#-------------------------------
rna_nascent <- read.table(opt$rna_nascent_fpkm, header=TRUE)
print("Determine expressed genes using 4sU data")
# what does the read count distribution look like?
# compute median per gene and plot it as histogram
rna_nascent$median_genecounts <- apply(rna_nascent[,-1], 1, FUN=median)
# lots of medians that are below 1
hist(log10(rna_nascent$median_genecounts))
# filter based on expression
expressed_genes <- rna_nascent %>% 
  dplyr::filter(median_genecounts > 0) 

#-------------------------------
### Load genekey to annotate ensembl to mgi
#-------------------------------
geneKey <- read.delim(opt$genekey)

#merge gene annotations to results table
expressed_genes <- merge(expressed_genes,
                         geneKey, 
                         by.x="Geneid", 
                         by.y="ensembl_gene_id")

#-------------------------------
## Getting the sequences
#-------------------------------
# import summit of ChIPseq peaks
ChIPseq_summits <- read.table(opt$chipseq_summits)
ChIPseq_ranges <- GRanges(seqnames = ChIPseq_summits[,c("V1")],
                          ranges = IRanges(start=ChIPseq_summits[,c("V2")],
                                           end=ChIPseq_summits[,c("V3")]-1)) # to make up for 0 vs 1 encoding
ChIPseq_ranges$id <- c(1:length(ChIPseq_ranges))


# NOTE: if we only want to annotate to genes that are expressed, we could use ChIPpeakAnno and a filtered annoDB object instead
# annotate it to genes
print("Annotate ChIPseq summit to closest gene (using genomic reference)")
summitAnno <- annotatePeak(ChIPseq_ranges, 
                           tssRegion=c(-3000, 3000), 
                           TxDb=txdb, annoDb = "org.Mm.eg.db")

summitAnno_df <- summitAnno %>% as.data.frame()

# see which ones are DE genes and add that info to GRanges as column "directionchange"
print("Add info on which genes are DE to the annotated summits")
DE_4sU <- read.delim(opt$contrast_DexVSDexLPS)

summitAnno_df <- left_join(summitAnno_df, 
                           DE_4sU[,c("mgi_symbol","padj","log2FoldChange")],
                           by = c("SYMBOL" = "mgi_symbol"))

summitAnno_df <- summitAnno_df %>% mutate(change = case_when(padj<0.05 & log2FoldChange > opt$log2fcthresh ~ "up",
                                                             padj<0.05& log2FoldChange < -opt$log2fcthresh ~"down",
                                                             TRUE ~ "ns")
)

# save info on gene annotation
ChIPseq_ranges$mgi_symbol[match(summitAnno_df$id, ChIPseq_ranges$id )] <- summitAnno_df$SYMBOL

# assign directionchange as metadata column
ChIPseq_ranges$directionchange[match(summitAnno_df$id, ChIPseq_ranges$id )] <- summitAnno_df$change

# ass distance to TSS
ChIPseq_ranges$distanceToTSS[match(summitAnno_df$id, ChIPseq_ranges$id )] <- summitAnno_df$distanceToTSS

#-------------------------------
## Prefilter motifdb to motifs that are expressed in celltype
#-------------------------------
print("Prefilter meme_db for those motifs expressed in our 4sU data")
meme_db <- read_meme(opt$meme_db_path) %>% 
  to_df()

meme_db_expressed <- meme_db %>% 
  # the altname slot of meme_db contains the gene symbol (this is database-specific)
  # avoid mismatches cased by casing and keep motif if at least one part of composite is expressed
  tidyr::separate(altname, into=c("tf1", "tf2"), sep="::",remove=FALSE) %>%
  filter( str_to_upper(tf1) %in% str_to_upper(expressed_genes$mgi_symbol) | str_to_upper(tf2) %in% str_to_upper(expressed_genes$mgi_symbol)) %>%
  # we don't need the split TF info downstream
  dplyr::select(!c("tf1","tf2"))


print("Number of motifs pre-filtering: ")
nrow(meme_db)
print("Number of motifs post-filter: ")
nrow(meme_db_expressed)

#-------------------------------
## OPTIONAL: only run with motifs of interest
#-------------------------------
meme_motifsOI <- 
  meme_db_expressed %>% 
  filter(
    grepl("STAT", str_to_upper(altname)) |
      grepl("NR3C", str_to_upper(altname))
  )

#-------------------------------
## FIGURES on peak gene annotation
#-------------------------------

# filter peaks for those annotated to genes that are expressed
summitAnno_expr <- subset(summitAnno,
                          summitAnno@anno$SYMBOL %in% expressed_genes$mgi_symbol)

# filter the df version in the same fashion
summitAnno_df_expr <- summitAnno_df %>% filter(SYMBOL %in% expressed_genes$mgi_symbol)

#---------------------------------
# --- some stats
#---------------------------------

distbygene <-  summitAnno_df_expr  %>% 
  group_by(SYMBOL, change) %>%
  summarise(min_dist=min(abs(distanceToTSS)), 
            mean_dist=mean(abs(distanceToTSS))) %>%
  ungroup() %>%
  mutate(logmindist=log2(min_dist+1))

# why 30kb cutoff
distbygene_allDE <- distbygene  %>% 
  filter(!change=="ns") %>%
  mutate(change=factor(change,levels=c("down","up")))

print("We need to justify why we picked a cutoff of 30kb.")
print("From a genecentric perspective, we want to include the peak regions that most likely have a regulating function on the gene.")

print("With a cutoff of 30kb, how many genes DONT have at least one peak within that range?")
tbl <- table(distbygene_allDE$min_dist > 30000)
tbl[2]/(tbl[1]+tbl[2])

print("How many genes do we lose of both sets by using that cutoff?")
print("In the upregulated fraction:")
table( (distbygene %>% filter(change=="up"))$min_dist > 30000)
print("In the downregulated fraction:")
table( (distbygene %>% filter(change=="down"))$min_dist > 30000)

print("Min and mean dist for the genes with log2FC >", opt$log2fcthresh)
distbygene %>%
  filter(change=="up") %>%
  summarise_all(mean) %>%
  print()
print("Min and mean dist for the genes with log2FC <", opt$log2fcthresh )
distbygene %>%
  filter(change=="down") %>%
  summarise_all(mean) %>%
  print()

print("How many peaks per UPregulated gene:")
summitAnno_df_expr  %>% 
  filter(change=="up")%>%
  filter(abs(distanceToTSS)<30000) %>%
  group_by(SYMBOL) %>%
  summarise(count=n()) %>%
  pull(count) %>% 
  mean()
print("How many peaks per DOWNregulated gene:")
summitAnno_df_expr  %>% 
  filter(change=="down")%>%
  filter(abs(distanceToTSS)<30000) %>%
  group_by(SYMBOL) %>%
  summarise(count=n()) %>%
  pull(count) %>% 
  mean()

print("The distances between the peaks mapping to the same gene.")
print("returns NA if only one 1 peak is annotated to the gene - those are excluded")
print("For upregulated genes:")
summitAnno_df_expr  %>% 
  filter(change=="up")%>%
  filter(abs(distanceToTSS)<30000) %>%
  group_by(SYMBOL) %>%
  summarise(meanpeakdist = mean(dist(distanceToTSS))) %>% # 
  filter(!is.na(meanpeakdist)) %>%
  pull(meanpeakdist) %>%
  mean()
print("For downregulated genes:")
summitAnno_df_expr  %>% 
  filter(change=="down")%>%
  filter(abs(distanceToTSS)<30000) %>%
  group_by(SYMBOL) %>%
  summarise(meanpeakdist = mean(dist(distanceToTSS))) %>% 
  filter(!is.na(meanpeakdist)) %>%
  pull(meanpeakdist) %>%
  mean()

#--------------------------------------
#- permutations
#--------------------------------------

# Difference in means
groupdiff <- diff(tapply(distbygene_allDE$min_dist, distbygene_allDE$change, mean))
print("The mean minimum distance is smaller for the upregulated set")
print(paste("The group difference is: ",groupdiff))
print("Permutation test to see if this difference between the groups is meaningful")

#Permutation test
permutation.test <- function(group, outcome, n, reference){
  distribution=c()
  result=0
  for(i in 1:n){
    distribution[i]=diff(by(outcome, sample(group, length(group), FALSE), mean))
  }
  result=sum(abs(distribution) >= abs(groupdiff))/(n)
  return(list(result, distribution, groupdiff))
}

permtest_res <- permutation.test(distbygene_allDE$change, distbygene_allDE$min_dist, 100000, groupdiff)

#--------------------------------------
#- export objects
#--------------------------------------

#---------------------------------
# --- export results of the permutation test
saveRDS(permtest_res,
        file=paste0(opt$outdir,"permtest_res.rds"))

#---------------------------------
# --- export up and downregulated summitfraction for deeptools

table(ChIPseq_ranges$directionchange)
dir.create(paste0(opt$outdir,"peaks_annot2DEgenes_30kb_log2FC0.58/"))
export.bed(ChIPseq_ranges %>% filter(directionchange == "up") %>% filter(abs(distanceToTSS)<30000) ,
           con=paste0(opt$outdir,"peaks_annot2DEgenes_30kb_log2FC0.58/UP_summit_unmerged.bed"))
export.bed(ChIPseq_ranges %>% filter(directionchange == "down") %>% filter(abs(distanceToTSS)<30000) ,
           con=paste0(opt$outdir,"peaks_annot2DEgenes_30kb_log2FC0.58/DOWN_summit_unmerged.bed"))

#-------------------------------
## export objects to run memes afterwards

saveRDS(ChIPseq_ranges,
        file=paste0(opt$outdir,"../memes_bioc/ChIPseq_summit_Granges.rds"))

saveRDS(meme_db_expressed,
        file=paste0(opt$outdir,"../memes_bioc/meme_db_4sUexpressed.rds"))

#-------------------------------
## export objects for ggplot figures
saveRDS(summitAnno,
        file=paste0(opt$outdir,"summitAnno.rds"))
saveRDS(summitAnno_expr,
        file=paste0(opt$outdir,"summitAnno_expr.rds"))
saveRDS(summitAnno_df_expr,
        file=paste0(opt$outdir,"summitAnno_df_expr.rds"))

sink()

