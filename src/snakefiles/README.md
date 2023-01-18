To keep the workflow organized, the individual datatypes are analyzed and integrated in several **snakemake subworkflows** located within `src/snakefiles/` that are included through the Snakefile in the project root directory.  

# Full list of subworkflows and rules 

## Snakefile (located within the project root)
* download_JASPAR
* make_bowtie_index
* get_reference_files
* get_remap_data
* gunzip_blacklist
* scan_nr3c1_genomewide

## Snakefile_utils.smk 
* generate_fastqc_report
* sort_raw_BAM_and_blacklist_filter
* filter_duplicates
* make_bw_from_bam
* make_forwardbw_from_bam
* make_reversebw_from_bam
* make_RPGC_bw_from_bam
* make_CPM_bw_from_bam
* generate_bed_from_idr
* sort_narrowPeak_bypval
* read_summit_from_IDR_peaks
* get_chromsize_from_fa
* get_chromsize_bed
* create_bed_interval_from_summit
* merge_overlapping_bed_interval
* bam_index

## Snakefile_ChIPseq.smk
* merge_GR2020
* symlink_GR2020
* symlink_histone
* trim_adapter_sequences
* multiqc
* fastq_align
* PE_fastq_align
* run_macs2_narrow
* run_macs2_broad
* run_macs2_GR2020
* GRChIP_run_idr_for_replicates
* merge_BAMs_from_replicates
* get_adjusted_libsize_atac
* get_adjusted_libsize_h3k27ac

## Snakefile_integrate_expression.smk 
* run_DE_analysis
* custom_libnorm
* make_atacdiff
* make_h3k27acdiff
* memes_summitregions
* memes_fimo_summitregions
* deeptools_computematrix_UPvsDOWN_chipseqonly
* deeptools_heatmap_UPvsDOWN_chipseqonly_png
* deeptools_heatmap_UPvsDOWN_chipseqonly_svg
* filter_fimo_for_motif
* deeptools_computematrix_at_motifofinterest
* deeptools_heatmap_at_motifofinterest

## Snakefile_abcmodel.smk
* download_BMDM_cage_data
* get_macrophagespecific_tss_from_cage
* download_juicer_tools
* prep_hiCdata
* compute_powerlaw_fit_from_hic
* sort_narrowPeak_file
* merge_BAM_from_ATACseq
* combine_tss_and_summitregions
* call_candidate_regions
* get_TPM_condition_expression
* quantify_enhancer_activity
* compute_ABC_scores
* filter_ABC_scores
* get_predictions
* run_fimo_enhancerregions

## Snakefile_figures.smk 
* visualize_DE_results
* figure_chipseq_prepdata
* figure_proxanno_prepdata
* figure_chipseq_plot
* figure_proxanno_plot
* abc_visualizations
* abc_predictions_prepdata
* abc_predictions
* figure_GLMs
* figure_supplemental_GLMs
* compute_tfactivity
* figure_tfactivity_etc
* quantify_features_and_label_distribution

