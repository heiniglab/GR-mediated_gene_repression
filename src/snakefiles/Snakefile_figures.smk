
rule visualize_DE_results:
    input:
        norm="results/current/rnaseq_4sU/DESeq_normalized_counts_nototal.tsv",
        metadata="results/current/rnaseq_4sU/anno_nototal.tsv",
        contrast="results/current/rnaseq_4sU/res_DexLPSvLPS_ext.tsv"
    output:      
        "results/current/Figures/Figure_rnaseq.pdf"
        # also creates other outputs (within params.outdir), but I don't bother listing them all individually
    params:
        outdir="results/current/rnaseq_4sU/figures",
        log2fcthresh = 0.58
    conda:
        "../conda/r41_env.yml"
    shell:
        """
        Rscript src/scripts/DE_visualizations_4sU.r -n {input.norm} -a {input.metadata} --contrast {input.contrast} \
        --log2fcthresh {params.log2fcthresh} -o {params.outdir}
        """

rule figure_chipseq_prepdata:
    input:
        bam="results/current/ChIP/GR_2020/bam/DexLPS_chipseq_GR_PE_merged_GRCm38_sorted_ftr_dedup.bam",
        summits="results/current/ChIP/GR_2020/idr/DexLPS_chipseq_GR_GRCm38_p0.1_idr0-05_summit",
        nr3c1fullsitematches="results/current/homer/nr3c1.fullsitematches.mm10.bed",
        nr3c1halfsitematches="results/current/homer/nr3c1.halfsitematches.mm10.bed",
    output:
        "results/current/integration/genomation/sm_summitranges.rds",
        "results/current/integration/genomation/sm_nr3c1fullsitematches_ranges.rds",
        "results/current/integration/genomation/sm_summitranges_w_nr3c1fullsitehit.rds",
        "results/current/integration/genomation/sm_nr3c1fullsitehits_within_summitranges.rds"
    params:
        outdir="results/current/integration/genomation/"
    conda:
        "../conda/r41_env.yml"
    shell:
        """
        Rscript src/scripts/figure_chipseq_prepdata.r --chipseq_bam {input.bam} --chipseq_summits {input.summits} \
        --nr3c1fullsitematches {input.nr3c1fullsitematches} --nr3c1halfsitematches {input.nr3c1halfsitematches} \
        -o {params.outdir}
        """

rule figure_proxanno_prepdata:
    input:
        summits="results/current/ChIP/GR_2020/idr/DexLPS_chipseq_GR_GRCm38_p0.1_idr0-05_summit",
        genekey="results/current/rnaseq_4sU/geneKey_biomart_mm_k100.txt",
        contrast="results/current/rnaseq_4sU/res_DexLPSvLPS_ext.tsv",
        meme_db="data/current/motifs/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt",
        rna_nascent_fpkm="data/current/rnaseq_4su/Count_matrix/FPKM.matrix.xls"
    output:
        # the first file here holds some logged numbers we might need in the manuscript
        "results/current/integration/figure_proxanno_prepdata.out",
        "results/current/integration/permtest_res.rds",
        "results/current/integration/peaks_annot2DEgenes_30kb_log2FC0.58/UP_summit_unmerged.bed",
        "results/current/integration/peaks_annot2DEgenes_30kb_log2FC0.58/DOWN_summit_unmerged.bed",
        "results/current/memes_bioc/ChIPseq_summit_Granges.rds",
        "results/current/memes_bioc/meme_db_4sUexpressed.rds",
        "results/current/integration/summitAnno.rds",
        "results/current/integration/summitAnno_expr.rds",
        "results/current/integration/summitAnno_df_expr.rds"
    params:
        log2fcthresh = 0.58,
        outdir="results/current/integration/"
    conda:
        "../conda/r41_env.yml"
    shell:
        """
        Rscript src/scripts/figure_proxanno_prepdata.r --log2fcthresh {params.log2fcthresh} --chipseq_summits {input.summits} --genekey {input.genekey} \
        --contrast_DexVSDexLPS {input.contrast} --meme_db_path {input.meme_db} --rna_nascent_fpkm {input.rna_nascent_fpkm} -o {params.outdir}
        """

rule figure_chipseq_plot:
    input:
        summitAnno="results/current/integration/summitAnno.rds",
        bed="results/current/ChIP/GR_2020/idr/DexLPS_chipseq_GR_GRCm38_p0.1_idr0-05.bed",
        summits="results/current/ChIP/GR_2020/idr/DexLPS_chipseq_GR_GRCm38_p0.1_idr0-05_summit",
        genomation_scorematrix="results/current/integration/genomation/sm_summitranges.rds",
        nr3c1fullsitematches="results/current/homer/nr3c1.fullsitematches.mm10.bed",
        nr3c1halfsitematches="results/current/homer/nr3c1.halfsitematches.mm10.bed",
        streme="results/current/memes_bioc/streme_100bp/streme.xml"
    output:      
        "results/current/Figures/Figure_chipseq.pdf"
    conda:
        "memes"
    shell:
        """
        Rscript src/scripts/figure_chipseq_plot.r --summitAnno {input.summitAnno} --chipseq_peaks {input.bed} --chipseq_summits {input.summits} \
        --nr3c1fullsitematches {input.nr3c1fullsitematches} --nr3c1halfsitematches {input.nr3c1halfsitematches} \
        --sm_summitranges {input.genomation_scorematrix} --streme {input.streme}
        """
#------------------
# incomplete
#------------------

rule figure_proxanno_plot:
    input:
        summitAnno_expr = "results/current/integration/summitAnno_expr.rds",
        summitAnno_df_expr = "results/current/integration/summitAnno_df_expr.rds",
        permtest_res = "results/current/integration/permtest_res.rds",
        fimo_results = "results/current/memes_bioc/fimo_1000bp/fimo.rds",
        chipseq_summit_granges = "results/current/memes_bioc/ChIPseq_summit_Granges.rds",
        deeptools = "results/current/integration/deeptools/GR_chipseq_UPandDOWN_libnorm_30kb_1000b_1000a_heatmap.png",
        deeptools_svg = "results/current/integration/deeptools/GR_chipseq_UPandDOWN_libnorm_30kb_1000b_1000a_heatmap.svg",
        streme_100bp_up = "results/current/memes_bioc/streme_100bp_up/streme.xml",
        streme_100bp_down = "results/current/memes_bioc/streme_100bp_down/streme.xml"
    output:
        "results/current/Figures/Figure_peakgeneannotation.pdf"
    conda:
        "memes"
    shell:
        """
        Rscript src/scripts/figure_proxanno_plot.r --summitAnno_expr {input.summitAnno_expr} --summitAnno_df_expr {input.summitAnno_df_expr} \
        --permtest_res {input.permtest_res} --fimo_results {input.fimo_results} --chipseq_summit_granges {input.chipseq_summit_granges} \
        --deeptools {input.deeptools} --streme_100bp_up {input.streme_100bp_up} --streme_100bp_down {input.streme_100bp_down}
        """

rule abc_visualizations:
    input:
        ABC_DexLPS_all = "results/current/abcmodel/abcscores_DexLPS_GRCm38/EnhancerPredictionsAllPutative_ftr0.02.tsv",
        ABC_LPS_all = "results/current/abcmodel/abcscores_LPS_GRCm38/EnhancerPredictionsAllPutative_ftr0.02.tsv",
        contrast_DexVSDexLPS = "results/current/rnaseq_4sU/res_DexLPSvLPS_ext.tsv",
        chipseq_summits = "results/current/ChIP/GR_2020/idr/DexLPS_chipseq_GR_GRCm38_p0.1_idr0-05_summit",
        igv = "results/current/abcmodel/igv_snapshot_cd83.png"
    output:
        "results/current/Figures/Figure_abcresults.pdf"
    conda:
        "../conda/r41_env.yml"
    shell:
        """
        Rscript src/scripts/abc_visualizations.r --ABC_DexLPS_all {input.ABC_DexLPS_all} --ABC_LPS_all {input.ABC_LPS_all} \
        --contrast_DexVSDexLPS {input.contrast_DexVSDexLPS} --chipseq_summits {input.chipseq_summits} --igv {input.igv}
        """

rule abc_predictions_prepdata:
    input:
        ABC_DexLPS_all = "results/current/abcmodel/abcscores_DexLPS_GRCm38/EnhancerPredictionsAllPutative_ftr0.02.tsv",
        ABC_LPS_all = "results/current/abcmodel/abcscores_LPS_GRCm38/EnhancerPredictionsAllPutative_ftr0.02.tsv",
        fimo_results_dexlps = "results/current/memes_bioc/fimo_DexLPSenhancers/fimo.rds",
        fimo_results_lps = "results/current/memes_bioc/fimo_LPSenhancers/fimo.rds",
        fimo_results_summitregion = "results/current/memes_bioc/fimo_1000bp/fimo.rds",
        chipseq_ranges = "results/current/memes_bioc/ChIPseq_summit_Granges.rds",
    output:
        # assignments
        "results/current/GLM_dataprep/assignment_summit_prox.rds",
        "results/current/GLM_dataprep/assignment_summits_ABCregion_dexlps.rds",
        "results/current/GLM_dataprep/assignment_summits_ABCregion_lps.rds",
        "results/current/GLM_dataprep/assignment_abcregion_dexlps.rds",
        "results/current/GLM_dataprep/assignment_abcregion_lps.rds",
        # motifcounts
        "results/current/GLM_dataprep/motifcounts_summitregion.rds",
        "results/current/GLM_dataprep/motifcounts_abcregion_dexlps.rds",
        "results/current/GLM_dataprep/motifcounts_abcregion_lps.rds",
    params: 
        outdir = "results/current/GLM_dataprep/"
    conda:
        "../conda/r41_env.yml"
    shell:
        """
        Rscript src/scripts/abc_predictions_prepdata.r --ABC_DexLPS_all {input.ABC_DexLPS_all} --ABC_LPS_all {input.ABC_LPS_all} \
        --fimo_results_dexlps {input.fimo_results_dexlps} --fimo_results_lps {input.fimo_results_lps} \
        --fimo_results_summitregion {input.fimo_results_summitregion} --chipseq_ranges {input.chipseq_ranges} --outdir {params.outdir}
        """

rule abc_predictions:
    input:
        contrast_DexLPSvLPS = "results/current/rnaseq_4sU/res_DexLPSvLPS_ext.tsv",
        # assignments
        assignment_summit_prox = "results/current/GLM_dataprep/assignment_summit_prox.rds",
        assignment_summits_abcregion_dexlps = "results/current/GLM_dataprep/assignment_summits_ABCregion_dexlps.rds",
        assignment_summits_abcregion_lps = "results/current/GLM_dataprep/assignment_summits_ABCregion_lps.rds",
        assignment_abcregion_dexlps = "results/current/GLM_dataprep/assignment_abcregion_dexlps.rds",
        assignment_abcregion_lps = "results/current/GLM_dataprep/assignment_abcregion_lps.rds",
        # motifcounts
        motifcounts_summitregion = "results/current/GLM_dataprep/motifcounts_summitregion.rds",
        motifcounts_abcregion_dexlps = "results/current/GLM_dataprep/motifcounts_abcregion_dexlps.rds",
        motifcounts_abcregion_lps = "results/current/GLM_dataprep/motifcounts_abcregion_lps.rds",
    output:
        "results/current/GLM_performance_scaled/gg_rawcounts.rds",
        "results/current/GLM_performance_scaled/model_coefs_joint.rds",
        "results/current/GLM_performance_scaled/model_coefs_sep.rds",
        "results/current/GLM_performance_scaled/AUC_metrics.rds"
    params:
        outdir="results/current/GLM_performance_scaled/"
    conda:
        "../conda/r41_env.yml"
    shell:
        """
        Rscript src/scripts/abc_predictions.r --contrast_DexLPSvLPS {input.contrast_DexLPSvLPS} \
        --assignment_summit_prox {input.assignment_summit_prox} --assignment_summits_abcregion_dexlps {input.assignment_summits_abcregion_dexlps} --assignment_summits_abcregion_lps {input.assignment_summits_abcregion_lps} \
        --assignment_abcregion_dexlps {input.assignment_abcregion_dexlps} --assignment_abcregion_lps {input.assignment_abcregion_lps} \
        --motifcounts_summitregion {input.motifcounts_summitregion} --motifcounts_abcregion_dexlps {input.motifcounts_abcregion_dexlps} --motifcounts_abcregion_lps {input.motifcounts_abcregion_lps} --outdir {params.outdir}
        """

rule figure_GLMs:
    input:
        model_coefs_joint = "results/current/GLM_performance_scaled/model_coefs_joint.rds",
        model_coefs_sep = "results/current/GLM_performance_scaled/model_coefs_sep.rds",
        auc = "results/current/GLM_performance_scaled/AUC_metrics.rds",
        motifcounts_summitregion = "results/current/GLM_dataprep/motifcounts_summitregion.rds",
        raw_counts = "results/current/GLM_performance_scaled/gg_rawcounts.rds"
    output:
        "results/current/Figures/Figure_motifcorrelations.pdf",
        "results/current/Figures/Figure_GLMs.pdf",
        "results/current/Figures/Figure_GLMs_coefssep.pdf"
    conda:
        "../conda/r41_env.yml"
    shell:
        """
        Rscript src/scripts/figure_GLMs.r --model_coefs_joint {input.model_coefs_joint} \
        --model_coefs_sep {input.model_coefs_sep} --auc {input.auc} \
        --motifcounts_summitregion {input.motifcounts_summitregion} --raw_counts {input.raw_counts}
        """

rule figure_supplemental_GLMs:
    input:
        auc_metrics = "results/current/GLM_performance_scaled/AUC_metrics.rds"
    params:
        dirname_featurematrizes = "results/current/GLM_featurematrizes_unscaled/",
        dirname_models = "results/current/GLM_performance_scaled/"
    output:
        png="results/current/Figures/Figure_supplemental_GLMs.png",
        pdf="results/current/Figures/Figure_supplemental_GLMs.pdf"
    conda:
        "../conda/r41_env.yml"
    shell:
        """
        Rscript src/scripts/figure_supplemental_GLMs.r --auc {input.auc_metrics} \
        --dirname_featurematrizes {params.dirname_featurematrizes} \
        --dirname_models {params.dirname_models} --outfig {output.png}
        """

rule compute_tfactivity:
    input:
         gene_annot = "data/current/gencode_annotations/gencode.vM23.annotation.gene.gtf",
         binding_sites_remap = "data/current/remap/remap2022_all_macs2_mm10_v1_0_separated_ftr.bed",
         rna_nascent_fpkm = "data/current/rnaseq_4su/Count_matrix/FPKM.matrix.xls",
         genekey = "results/current/rnaseq_4sU/geneKey_biomart_mm_k100.txt",
         cage = "results/current/abcmodel/tss_and_genecoords/mac_cage_maxscore.bed"
    output:
         "results/current/tfactivity/TFA_4sU.tsv",
         "results/current/tfactivity/FPKMcounts_mgiaggr.tsv",
         "results/current/tfactivity/TFA_4sU_CAGE.tsv"
    conda:
        "../conda/r41_env.yml"
    params:
        outdir = "results/current/tfactivity/"
    shell:
        """
        Rscript src/scripts/tf_activity_4sU.r --gene_annot {input.gene_annot} \
        --binding_sites_remap {input.binding_sites_remap} --rna_nascent_fpkm {input.rna_nascent_fpkm} \
        --genekey {input.genekey} --cage {input.cage} --outdir {params.outdir}
        """

rule figure_tfactivity_etc:
    input:
         # footprint profileplots and IGV screenshots are currently added in manually
         tfactivity = "results/current/tfactivity/TFA_4sU_CAGE.tsv",
         expr = "results/current/tfactivity/FPKMcounts_mgiaggr.tsv",
         difffootprint = "data/current/atacseq/footprints/DexLPSvLPS_diff_footprint/differential_statistics.txt",
         memedb_expressed = "results/current/memes_bioc/meme_db_4sUexpressed.rds",
         heatmap = "results/current/integration/deeptools/bw_at_STAT3_heatmap_cpm.png",
         chipms = "data/current/chipMS_macro/KopievonGRChIPMS.xlsx"
    output:
        "results/current/Figures/Figure_stats.pdf",
        "results/current/Figures/Suppl_Figure_statmotifs.pdf"
    conda:
        "../conda/r41_env.yml"
    shell:
        """
        Rscript src/scripts/figure_stats.r --tfactivity {input.tfactivity} --expr {input.expr} --difffootprint {input.difffootprint} \
        --memedb_expressed {input.memedb_expressed} --heatmap {input.heatmap} --chipms {input.chipms}
        """

rule quantify_features_and_label_distribution:
    input:
        contrast_DexLPSvLPS = "results/current/rnaseq_4sU/res_DexLPSvLPS_ext.tsv",
        # assignments
        assignment_summit_prox = "results/current/GLM_dataprep/assignment_summit_prox.rds",
        assignment_summits_abcregion_dexlps = "results/current/GLM_dataprep/assignment_summits_ABCregion_dexlps.rds",
        assignment_summits_abcregion_lps = "results/current/GLM_dataprep/assignment_summits_ABCregion_lps.rds",
        assignment_abcregion_dexlps = "results/current/GLM_dataprep/assignment_abcregion_dexlps.rds",
        assignment_abcregion_lps = "results/current/GLM_dataprep/assignment_abcregion_lps.rds",
        # motifcounts
        motifcounts_summitregion = "results/current/GLM_dataprep/motifcounts_summitregion.rds",
        motifcounts_abcregion_dexlps = "results/current/GLM_dataprep/motifcounts_abcregion_dexlps.rds",
        motifcounts_abcregion_lps = "results/current/GLM_dataprep/motifcounts_abcregion_lps.rds",
        # coefficients of selected models
        model_coefs_joint = "results/current/GLM_performance_scaled/model_coefs_joint.rds",
        model_coefs_sep = "results/current/GLM_performance_scaled/model_coefs_sep.rds"
    output:
        "results/current/Figures/metrics_n_features_label_distro.tsv"
    params:
        featuredir = "results/current/GLM_featurematrizes_unscaled/"
    conda:
        "../conda/r41_env.yml"
    shell:
        """
        Rscript src/scripts/quantify_features_and_label_distribution.r --contrast_DexLPSvLPS {input.contrast_DexLPSvLPS} \
        --assignment_summit_prox {input.assignment_summit_prox} \
        --assignment_summits_abcregion_dexlps {input.assignment_summits_abcregion_dexlps} \
        --assignment_summits_abcregion_lps {input.assignment_summits_abcregion_lps} \
        --assignment_abcregion_dexlps {input.assignment_abcregion_dexlps} \
        --assignment_abcregion_lps {input.assignment_abcregion_lps} \
        --motifcounts_summitregion {input.motifcounts_summitregion} \
        --motifcounts_abcregion_dexlps {input.motifcounts_abcregion_dexlps} \
        --motifcounts_abcregion_lps {input.motifcounts_abcregion_lps} \
        --model_coefs_joint {input.model_coefs_joint} \
        --model_coefs_sep {input.model_coefs_sep} \
        --featuredir {params.featuredir} --outfile {output}
        """  
