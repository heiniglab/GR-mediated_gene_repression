###############################################################################
#---------------------define snakefiles that create prerequisites
###############################################################################

include: "Snakefile_utils.smk"
        
###############################################################################
#---------------------run downstream analyses
###############################################################################

# This workflow depends on the output of the bioinformatics RNAseq pipeline, cloned from `https://ascgitlab.helmholtz-muenchen.de/bioinformatics_core/pipeline_RNA-Seq.git`

rule run_DE_analysis:
    input:
        counts="data/current/rnaseq_4su/Count_matrix/Count.matrix.xls",
        metadata="data/current/rnaseq_4su/raw_fastq/metadata_4sU.txt"
    output:
        # normalized counts
        "results/current/rnaseq_4sU/DESeq_normalized_counts_nototal.tsv",
        "results/current/rnaseq_4sU/DESeq_normalized_counts_nototal.rds",
        # properly formatted sample annotations
        "results/current/rnaseq_4sU/anno_nototal.tsv",
        # pairwise contrats
        "results/current/rnaseq_4sU/contrast_DexVSDexLPS.rds",
        "results/current/rnaseq_4sU/contrast_LPSvVeh.rds",
        "results/current/rnaseq_4sU/contrast_DexLPSvVeh.rds",
        # the biomart version used to annotate ensembl IDs to MGI symbols
        "results/current/rnaseq_4sU/geneKey_biomart_mm_k100.txt",
        # results as table with MGI annotation
        "results/current/rnaseq_4sU/res_DexLPSvLPS_ext.tsv",
        "results/current/rnaseq_4sU/res_LPSvVeh_ext.tsv",
        "results/current/rnaseq_4sU/res_DexLPSvVeh_ext.tsv"
    params:
        outdir="results/current/rnaseq_4sU"
    conda:
        "../conda/r41_env.yml"
    shell:
        """
        mkdir -p {params.outdir} && \
        Rscript src/scripts/DE_analysis_4sU.r -c {input.counts} -a {input.metadata} -k 100 -o {params.outdir}
        """

rule custom_libnorm:
    input:
        bw_DexLPS_h3k27ac = "results/current/ChIP/H3K27ac/bw/DexLPS_histone_H3K27ac_PE_merged_GRCm38.bw",
        bw_LPS_h3k27ac = "results/current/ChIP/H3K27ac/bw/LPS_histone_H3K27ac_PE_merged_GRCm38.bw",
        counts_h3k27ac = "results/current/ChIP/H3K27ac/counts",
        bw_DexLPS_atac = "results/current/atacseq/bw/merged_DexLPS_GRCm38.bw",
        bw_LPS_atac = "results/current/atacseq/bw/merged_LPS_GRCm38.bw",
        counts_atac = "results/current/atacseq/counts",
        gtf = "data/current/Mus_musculus.GRCm38.100.gtf.gz"
    output:
        "results/current/ChIP/H3K27ac/bw/DexLPS_histone_H3K27ac_PE_merged_GRCm38_libnorm.bw",
        "results/current/ChIP/H3K27ac/bw/LPS_histone_H3K27ac_PE_merged_GRCm38_libnorm.bw",
        "results/current/atacseq/bw/merged_DexLPS_GRCm38_libnorm.bw",
        "results/current/atacseq/bw/merged_LPS_GRCm38_libnorm.bw",
        "data/current/Mus_musculus.GRCm38.100_ftrnogene.gtf",
        "data/current/Mus_musculus.GRCm38.100_ftrnogene.gff3"
    conda:
        "../conda/r41_env.yml"
    shell:
        """
        Rscript src/scripts/custom_bw_libnorm.r \
        --bw_DexLPS_h3k27ac {input.bw_DexLPS_h3k27ac} --bw_LPS_h3k27ac {input.bw_LPS_h3k27ac} --counts_h3k27ac {input.counts_h3k27ac} \
        --bw_DexLPS_atac {input.bw_DexLPS_atac} --bw_LPS_atac {input.bw_LPS_atac} --counts_atac {input.counts_atac} --gtf {input.gtf}
        """

rule make_atacdiff:
    input:
        bw1="results/current/atacseq/bw/merged_DexLPS_GRCm38_libnorm.bw",
        bw2="results/current/atacseq/bw/merged_LPS_GRCm38_libnorm.bw"
    output:
        "results/current/atacseq/bw/merged_GRCm38_libnorm_diff.bw"
    conda:
        "../conda/deeptools.yml"
    shell:
        """
        bigwigCompare --bigwig1 {input.bw1} --bigwig2 {input.bw2} -p 4 -o {output}
        """

rule make_h3k27acdiff:
    input:
        bw1="results/current/ChIP/H3K27ac/bw/DexLPS_histone_H3K27ac_PE_merged_GRCm38_libnorm.bw",
        bw2="results/current/ChIP/H3K27ac/bw/LPS_histone_H3K27ac_PE_merged_GRCm38_libnorm.bw",
    output:
        "results/current/ChIP/H3K27ac/bw/histone_H3K27ac_PE_merged_GRCm38_libnorm_diff.bw",
    conda:
        "../conda/deeptools.yml"
    shell:
        """
        bigwigCompare --bigwig1 {input.bw1} --bigwig2 {input.bw2} -p 4 -o {output}
        """

# the streme analyses are time intensive and only get run, if the output file does not exist already
rule memes_summitregions:
    input:
        summit_granges="results/current/memes_bioc/ChIPseq_summit_Granges.rds",
        memedb_expressed="results/current/memes_bioc/meme_db_4sUexpressed.rds",
    output:
        # separate for up- and down and smaller region
        "results/current/memes_bioc/streme_100bp_down/streme.xml",
        "results/current/memes_bioc/streme_100bp_up/streme.xml",
        "results/current/memes_bioc/dreme_100bp_down/dreme.xml",
        "results/current/memes_bioc/dreme_100bp_up/dreme.xml",
        # discriminative:
        # enriched in up with down as control
        "results/current/memes_bioc/ame_discr_up/ame.tsv",
        # enriched in down with up as control
        "results/current/memes_bioc/ame_discr_down/ame.tsv",
        # all and big region
        #"results/current/memes_bioc/streme_1000bp/streme.xml",
        done="results/current/memes_bioc/done"
    conda:
        "memes"
    shell:
        """
        mkdir -p results/current/memes_bioc/ && \
        perlbrew use perl-5.34.0 && \
        Rscript src/scripts/memes_runanalyses.r \
        --summit_granges {input.summit_granges} --memedb_expressed {input.memedb_expressed} && \
        touch {output.done}
        """

# the streme analyses are time intensive and only get run, if the output file does not exist already
rule memes_fimo_summitregions:
    input:
        summit_granges="results/current/memes_bioc/ChIPseq_summit_Granges.rds",
        memedb_expressed="results/current/memes_bioc/meme_db_4sUexpressed.rds",
    output:
        # fimo hits of motifs from expressed TFs
        "results/current/memes_bioc/fimo_1000bp/fimo.rds"
    params:
    conda:
        "../conda/r41_env.yml"
    shell:
        """
        mkdir -p $(dirname {output}) && \
        perlbrew use perl-5.34.0 && \
        Rscript src/scripts/memes_fimo_runanalyses.r \
        --summit_granges {input.summit_granges} --memedb_expressed {input.memedb_expressed}
        """

rule deeptools_computematrix_UPvsDOWN_chipseqonly:
    input:
        tUP="results/current/integration/peaks_annot2DEgenes_30kb_log2FC0.58/UP_summit_unmerged.bed",
        tDOWN="results/current/integration/peaks_annot2DEgenes_30kb_log2FC0.58/DOWN_summit_unmerged.bed",
        chipseq="results/current/ChIP/GR_2020/bw/DexLPS_chipseq_GR_PE_merged_GRCm38_cpm.bw",
        atac_diff="results/current/atacseq/bw/merged_GRCm38_libnorm_diff.bw",
        h3k27ac_diff="results/current/ChIP/H3K27ac/bw/histone_H3K27ac_PE_merged_GRCm38_libnorm_diff.bw"
    output:
        "results/current/integration/deeptools/GR_chipseq_UPandDOWN_libnorm_30kb_1000b_1000a_matrix.gz"
    conda:
        "../conda/deeptools.yml"
    shell:
        """
        mkdir -p "results/current/integration/deeptools" &&
        computeMatrix reference-point --referencePoint center -b 1000 -a 1000 -bs 50 -R {input.tUP} {input.tDOWN} \
        -S {input.chipseq} {input.atac_diff} {input.h3k27ac_diff} -o {output}
        """

rule deeptools_heatmap_UPvsDOWN_chipseqonly_png:
    input:
        "results/current/integration/deeptools/GR_chipseq_UPandDOWN_libnorm_30kb_1000b_1000a_matrix.gz"
    output:
        "results/current/integration/deeptools/GR_chipseq_UPandDOWN_libnorm_30kb_1000b_1000a_heatmap.png"
    conda:
        "../conda/deeptools.yml"
    shell:
        """
        plotHeatmap -m {input} \
        --sortRegions descend --refPointLabel summit \
        --samplesLabel chipseq_DexLPS atac_diff h3k27ac_diff \
        --colorMap Greens RdYlBu_r RdYlBu_r \
        --zMin 0 -2 -2  --zMax 0.8 1.8 1.8 \
        --regionsLabel "regions of upregulated genes" "regions of downregulated genes" \
        -out {output}
        """

rule deeptools_heatmap_UPvsDOWN_chipseqonly_svg:
    input:
        "results/current/integration/deeptools/GR_chipseq_UPandDOWN_libnorm_30kb_1000b_1000a_matrix.gz"
    output:
        "results/current/integration/deeptools/GR_chipseq_UPandDOWN_libnorm_30kb_1000b_1000a_heatmap.svg"
    conda:
        "../conda/deeptools.yml"
    shell:
        """
        plotHeatmap -m {input} \
        --sortRegions descend --refPointLabel summit \
        --samplesLabel chipseq_DexLPS atac_diff h3k27ac_diff \
        --colorMap Greens RdYlBu_r RdYlBu_r \
        --zMin 0 -2 -2  --zMax 0.8 1.8 1.8 \
        --regionsLabel "regions of upregulated genes" "regions of downregulated genes" \
        -out {output}
        """

rule filter_fimo_for_motif:
    input:
        fimo = "results/current/memes_bioc/fimo_1000bp/fimo.rds",
        summit_granges = "results/current/memes_bioc/ChIPseq_summit_Granges.rds"
    output:
        outfile = "results/current/memes_bioc/{motif_altname}_hits_100bp.bed"
    conda:
        "../conda/r41_env.yml"
    shell:
        """
        Rscript src/scripts/filter_fimo_for_motif.r --fimo {input.fimo} --summit_granges {input.summit_granges} \
        --motif_altname {wildcards.motif_altname}
        """

rule deeptools_computematrix_at_motifofinterest:
    input:
        #created with "filter_fimo_for_motif.r"
        targets="results/current/memes_bioc/STAT3_hits_100bp.bed",
        chipseq="results/current/ChIP/GR_2020/bw/DexLPS_chipseq_GR_PE_merged_GRCm38_cpm.bw",
        atac_dexlps="results/current/atacseq/bw/merged_DexLPS_GRCm38_libnorm.bw",
        atac_lps="results/current/atacseq/bw/merged_LPS_GRCm38_libnorm.bw",
        h3k27ac_dexlps="results/current/ChIP/H3K27ac/bw/DexLPS_histone_H3K27ac_PE_merged_GRCm38_libnorm.bw",
        h3k27ac_lps="results/current/ChIP/H3K27ac/bw/LPS_histone_H3K27ac_PE_merged_GRCm38_libnorm.bw",
    output:
        "results/current/integration/deeptools/bw_at_STAT3_matrix_cpm.gz", 
    conda:
        "../conda/deeptools.yml"
    shell:
        """
        mkdir -p $(dirname {output}) && \
        computeMatrix reference-point --referencePoint center -b 500 -a 500 -bs 2 -R {input.targets} \
        -S {input.chipseq} {input.atac_dexlps} {input.atac_lps} {input.h3k27ac_dexlps} {input.h3k27ac_lps} -o {output}
        """

rule deeptools_heatmap_at_motifofinterest:
    input:
        "results/current/integration/deeptools/bw_at_STAT3_matrix_cpm.gz", 
    output:
        "results/current/integration/deeptools/bw_at_STAT3_heatmap_cpm.png", 
    conda:
        "../conda/deeptools.yml"
    shell:
        """
        plotHeatmap -m {input} \
        --sortRegions descend --refPointLabel summit \
        --samplesLabel chipseq_DexLPS atac_dexlps atac_lps h3k27ac_dexlps h3k27ac_lps \
        --colorMap Greens RdYlBu_r RdYlBu_r RdYlBu_r RdYlBu_r \
        --zMin 0 0 0 0 0  --zMax 0.8 10 10 10 10  \
        -out {output}
        """
