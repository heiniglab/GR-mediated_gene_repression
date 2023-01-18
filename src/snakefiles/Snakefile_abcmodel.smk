
# git clone https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction.git retreived on July 6th 2021

include: "Snakefile_utils.smk"

###############################################################################
#---------------------Step 1. Define candidate elemets
###############################################################################

rule download_BMDM_cage_data:
    input:
    output:
        s1 = "data/current/FANTOM5_CAGE/BMDM_pool1.ctss",
        s2 = "data/current/FANTOM5_CAGE/BMDM_pool2.ctss"
    params:
        link1 = "http://fantom.gsc.riken.jp/5/datafiles/latest/basic/mouse.timecourse.hCAGE/macrophage%252c%2520bone%2520marrow%2520derived%252c%2520pool1.CNhs11457.3560-170A1.mm9.ctss.bed.gz",
        link2 = "http://fantom.gsc.riken.jp/5/datafiles/latest/basic/mouse.timecourse.hCAGE/macrophage%252c%2520bone%2520marrow%2520derived%252c%2520pool2.CNhs11532.3632-171A1.mm9.ctss.bed.gz",
        temp_s1 = "data/current/FANTOM5_CAGE/macrophage\,\ bone\ marrow\ derived\,\ pool1.CNhs11457.3560-170A1.mm9.ctss.bed.gz", #macrophage, bone marrow derived, pool1.CNhs11457.3560-170A1.mm9.ctss.bed.gz",
        temp_s2 = "data/current/FANTOM5_CAGE/macrophage\,\ bone\ marrow\ derived\,\ pool2.CNhs11532.3632-171A1.mm9.ctss.bed.gz", #macrophage, bone marrow derived, pool2.CNhs11532.3632-171A1.mm9.ctss.bed.gz",
        outdir = "data/current/FANTOM5_CAGE/"
    shell:
        """
        mkdir -p {params.outdir} && \
        wget {params.link1} -P {params.outdir} && \
        zcat {params.temp_s1} | awk '{{print $1,$2,$6,$5}}' OFS='\t' - > {output.s1} && \
        wget {params.link2} -P {params.outdir} && \
        zcat {params.temp_s2} | awk '{{print $1,$2,$6,$5}}' OFS='\t' - > {output.s2}
        """

#TSS gets used in "includeregions" argument, genecoords in Step 2: quantify enhancer activity
rule get_macrophagespecific_tss_from_cage:
    input:
        #"results/current/setup_cager_env.done",
        ctss_pool1 = "data/current/FANTOM5_CAGE/BMDM_pool1.ctss",
        ctss_pool2 = "data/current/FANTOM5_CAGE/BMDM_pool2.ctss",
        liftoverchain = "data/current/mm9ToMm10.over.chain", 
        gencode_mm9_geneanno = "data/current/gencode_annotations/gencode.vM1.annotation.gtf.gz",
        gencode_mm10_geneanno = "data/current/gencode_annotations/gencode.vM23.annotation.gene.gtf"
    output:
        "results/current/abcmodel/tss_and_genecoords/mac_cage_tssclusterregions.bed",
        "results/current/abcmodel/tss_and_genecoords/mac_cage_dominant_ctss.bed",
        "results/current/abcmodel/tss_and_genecoords/mac_cage_maxscore.bed",
        "results/current/abcmodel/tss_and_genecoords/reference_genecoords.bed",
        "results/current/abcmodel/tss_and_genecoords/reference_promoterregions.bed"
    params:
        outdir="results/current/abcmodel/tss_and_genecoords/"
    conda:
        #"cager" use in case local env should be used (e.g. due to postlink script issues of go.db and org.mm.eg.db)
        "../conda/cager.yml" 
    shell:
        """
        mkdir -p {params.outdir} && \
        Rscript src/scripts/get_TSS_from_cage.r --ctss_pool1 {input.ctss_pool1} --ctss_pool2 {input.ctss_pool2} \
        --liftoverchain {input.liftoverchain} --gencode_mm9_geneanno {input.gencode_mm9_geneanno} \
        --gencode_mm10_geneanno {input.gencode_mm10_geneanno} --outdir {params.outdir}
        """

rule download_juicer_tools:
    input:
    output:
      "src/juicer_tools.1.9.9_jcuda.0.8.jar"
    conda:
        "../conda/abcenv.yml"      
    shell:
        """
        wget https://hicfiles.tc4ga.com/public/juicer/juicer_tools.1.9.9_jcuda.0.8.jar -P src/ &&
        ln -s juicer_tools.1.9.9_jcuda.0.8.jar  src/juicer_tools.jar
        """ 

rule prep_hiCdata:
    input:
      hic="data/current/hic/GSM3181429_BMDM_HiC_B1.hic",
      jar="src/juicer_tools.jar"
    output:
      "results/current/abcmodel/hic/raw/chr1/chr1.KRnorm.gz" # TODO: rewrite to include all chromosomes?
    conda:
      "../conda/abcenv.yml"
    params:
      outdir="results/current/abcmodel/hic/raw/"
    shell:
      """
      mkdir -p {params.outdir} && \
      python src/scripts/abcmodel/juicebox_dump.py \
      --hic_file {input.hic} \
      --juicebox "java -jar src/juicer_tools.jar" \
      --outdir {params.outdir} \
      --chromosomes 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X
      """

rule compute_powerlaw_fit_from_hic:
  input:
    "results/current/abcmodel/hic/raw/chr1/chr1.KRnorm.gz"
  output:
    "results/current/abcmodel/hic/raw/powerlaw/hic.powerlaw.txt"
  conda:
    "../conda/abcenv.yml"
  params:
    hicdir="results/current/abcmodel/hic/raw/",
    outdir="results/current/abcmodel/hic/raw/powerlaw/"
  shell:
    """
    python src/scripts/abcmodel/compute_powerlaw_fit_from_hic.py \
    --hicDir {params.hicdir} \
    --outDir {params.outdir} \
    --maxWindow 1000000 \
    --minWindow 5000 \
    --resolution 5000 \
    --chr $(echo chr{{1..19}} chrX | tr ' ' ,)
    """
 
rule sort_narrowPeak_file:
    input:
        chromsizes="results/current/genomesize/GRCm38.chromsizes",
        narrowpeak="data/current/atacseq/merged/{condition}/peaks.bed"
    conda:
        "../conda/samtools_19.yml"
    output:
        "results/current/atacseq/macs2/{condition}_peaks.narrowPeak.sorted" 
    shell:
        """
        mkdir -p results/current/atacseq/macs2/ && \
        bedtools sort -faidx {input.chromsizes} -i {input.narrowpeak} > {output}
        """

rule merge_BAM_from_ATACseq:
    input:
        "data/current/atacseq/bam/processed/{condition}", # this is a directory
    output:
        "results/current/atacseq/bam/merged_{condition}_GRCm38_algnd.bam", # ending in _algnd.bam lets it get recognized by the general rule
    conda:
        "../conda/samtools_19.yml"
    shell:
        "samtools merge {output} {input}/*.bam"

rule combine_tss_and_summitregions:
    input:
        mac_cage_tssclusterregions = "results/current/abcmodel/tss_and_genecoords/mac_cage_tssclusterregions.bed",
        ref_promoterregions = "results/current/abcmodel/tss_and_genecoords/reference_promoterregions.bed",
        summitintervals_GR_ChIPseq = "results/current/ChIP/GR_2020/idr/DexLPS_chipseq_GR_GRCm38_p0.1_idr0-05_summitinterval_slop250_merged.bed"
    output:
        "results/current/abcmodel/includeregions.bed"
    shell:
        "cat {input.mac_cage_tssclusterregions} {input.ref_promoterregions} {input.summitintervals_GR_ChIPseq} > {output}"

rule call_candidate_regions:
    input:
        narrowpeak_sorted="results/current/atacseq/macs2/{condition}_peaks.narrowPeak.sorted", 
        atac_bam="results/current/atacseq/bam/merged_{condition}_GRCm38_sorted_ftr_dedup.bam",
        atac_bai="results/current/atacseq/bam/merged_{condition}_GRCm38_sorted_ftr_dedup.bai",
        chromsizes="results/current/genomesize/GRCm38.chromsizes",
        chromsizes_bed="results/current/genomesize/GRCm38.chromsizes.bed",
        blacklist="data/current/mm10.blacklist.bed",
        includeregions="results/current/abcmodel/includeregions.bed"
    conda:
        "../conda/abcenv.yml"
    params:
        outdir="results/current/abcmodel/candidateregions/"
    output:
        "results/current/abcmodel/candidateregions/{condition}_peaks.narrowPeak.sorted.candidateRegions.bed" 
        # filtered, extended and merged peak calls from MACS2. These are the candidate regions used in downstream scripts.
    shell:
        """
        python src/scripts/abcmodel/makeCandidateRegions.py \
        --narrowPeak {input.narrowpeak_sorted} \
        --bam {input.atac_bam} \
        --outDir {params.outdir} \
        --chrom_sizes {input.chromsizes} \
        --regions_blocklist {input.blacklist} \
        --regions_includelist {input.includeregions} \
        --nStrongestPeaks 150000
        """
# 
# ###############################################################################
# #---------------------Step 2. Quantifying Enhancer Activity
# ###############################################################################
# 

rule get_TPM_condition_expression:
    input:
        tpm="data/current/rnaseq_4su/Count_matrix/TPM.matrix.xls",
        genekey="results/current/rnaseq_4sU/geneKey_biomart_mm_k100.txt"
    output:
        "results/current/abcmodel/expression/Veh_tpm.tsv",
        "results/current/abcmodel/expression/DexLPS_tpm.tsv",
        "results/current/abcmodel/expression/LPS_tpm.tsv"
    conda:
        "../conda/r41_env.yml"
    shell:
        """
        mkdir -p "results/current/abcmodel/expression/" && \
        Rscript src/scripts/get_TPM_condition_expression.r -g {input.genekey} -e {input.tpm}
        """

rule quantify_enhancer_activity:
    input:
        candidateRegions="results/current/abcmodel/candidateregions/{condition}_peaks.narrowPeak.sorted.candidateRegions.bed" ,
        genecoords="results/current/abcmodel/tss_and_genecoords/reference_genecoords.bed",
        H3K27ac_bam1="results/current/ChIP/H3K27ac/bam/{condition}_histone_H3K27ac_rep1_PE_{basename}_sorted_ftr_dedup.bam",
        H3K27ac_bam2="results/current/ChIP/H3K27ac/bam/{condition}_histone_H3K27ac_rep2_PE_{basename}_sorted_ftr_dedup.bam",
        H3K27ac_bai1="results/current/ChIP/H3K27ac/bam/{condition}_histone_H3K27ac_rep1_PE_{basename}_sorted_ftr_dedup.bai",
        H3K27ac_bai2="results/current/ChIP/H3K27ac/bam/{condition}_histone_H3K27ac_rep2_PE_{basename}_sorted_ftr_dedup.bai",
        atac_bam="results/current/atacseq/bam/merged_{condition}_GRCm38_sorted_ftr_dedup.bam",
        chromsizes="results/current/genomesize/GRCm38.chromsizes",
        expression_tpm="results/current/abcmodel/expression/{condition}_tpm.tsv"
#         UbiquitouslyExpressedGenes="",
    output:
        enhancerlist="results/current/abcmodel/enhanceractivity_{condition}_{basename}/EnhancerList.txt", # Candidate enhancer regions with Dnase-seq and H3K27ac ChIP-seq read counts
        enhancerlist_bed="results/current/abcmodel/enhanceractivity_{condition}_{basename}/EnhancerList.bed",
        genelist="results/current/abcmodel/enhanceractivity_{condition}_{basename}/GeneList.txt", # Dnase-seq and H3K27ac ChIP-seq read counts on gene bodies and gene promoter regions
    conda:
        "../conda/abcenv.yml"
    params:
        outdir="results/current/abcmodel/enhanceractivity_{condition}_{basename}/",
        celltype="macrophage"
    shell:
#        --ubiquitously_expressed_genes {input.UbiquitouslyExpressedGenes} \
        """
        python src/scripts/abcmodel/run.neighborhoods.py \
        --candidate_enhancer_regions {input.candidateRegions} \
        --genes {input.genecoords} \
        --expression_table {input.expression_tpm} \
        --H3K27ac {input.H3K27ac_bam1},{input.H3K27ac_bam1} \
        --DHS {input.atac_bam} \
        --chrom_sizes {input.chromsizes} \
        --cellType {params.celltype} \
        --outdir {params.outdir} 
        """
# 
# ###############################################################################
# #---------------------Step 3. Computing the ABC Score
# ###############################################################################
# 
# # Compute ABC scores by combining Activity (as calculated by run.neighborhoods.py) and Hi-C.
#  
rule compute_ABC_scores:
    input:
        enhancerlist="results/current/abcmodel/enhanceractivity_{condition}_{basename}/EnhancerList.txt",
        genelist="results/current/abcmodel/enhanceractivity_{condition}_{basename}/GeneList.txt",
        chromsizes="results/current/genomesize/{basename}.chromsizes",
        hic="results/current/abcmodel/hic/raw/powerlaw/hic.powerlaw.txt"
    output:
       "results/current/abcmodel/abcscores_{condition}_{basename}/EnhancerPredictionsAllPutative.txt.gz"
    conda:
        "../conda/abcenv.yml"
    params:
        hic_dir="results/current/abcmodel/hic/raw",
        outdir="results/current/abcmodel/abcscores_{condition}_{basename}",
        celltype="macrophage",
        hic_res=5000
    shell:
        """
        python src/scripts/abcmodel/predict.py \
        --enhancers {input.enhancerlist} \
        --genes {input.genelist} \
        --HiCdir {params.hic_dir} \
        --chrom_sizes {input.chromsizes} \
        --hic_resolution {params.hic_res} \
        --scale_hic_using_powerlaw \
        --threshold .02 \
        --cellType {params.celltype} \
        --outdir {params.outdir} \
        --make_all_putative \
        --chromosomes $(echo chr{{1..19}} chrX | tr ' ' ,)
        """

rule filter_ABC_scores:
    input:
        "results/current/abcmodel/abcscores_{condition}_{basename}/EnhancerPredictionsAllPutative.txt.gz"
    output:
       "results/current/abcmodel/abcscores_{condition}_{basename}/EnhancerPredictionsAllPutative_ftr0.02.tsv"
    conda:
        "../conda/abcenv.yml"
    shell:
        """
        zcat {input} |awk '{{if ($21>=0.02 && $21!="NaN") {{print}}}}' FS='\\t' OFS='\t' > {output}
        """

###############################################################################
#--------------Make featurematrix using active regions as input regions
###############################################################################


rule run_fimo_enhancerregions:
    input:
        abc = "results/current/abcmodel/abcscores_{condition}_GRCm38/EnhancerPredictionsAllPutative_ftr0.02.tsv",
        memedb = "results/current/memes_bioc/meme_db_4sUexpressed.rds",
    output:
        rds_fimo = "results/current/memes_bioc/fimo_{condition}enhancers/fimo.rds",
    conda:
        "../conda/r41_env.yml"
    shell:
        """
        mkdir -p $(dirname {output}) && \
        perlbrew use perl-5.34.0 && \
        Rscript src/scripts/memes_runanalyses_ABCenhancerregions.r \
        --ABC_all {input.abc} --memedb_expressed {input.memedb} --output {output.rds_fimo}
        """

