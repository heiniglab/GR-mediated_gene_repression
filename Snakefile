
###############################################################################
#-------------------define snakefiles that create prerequisites
###############################################################################

include: "src/snakefiles/Snakefile_abcmodel.smk"
include: "src/snakefiles/Snakefile_ChIPseq.smk"
include: "src/snakefiles/Snakefile_figures.smk"
include: "src/snakefiles/Snakefile_integrate_expression.smk"
include: "src/snakefiles/Snakefile_utils.smk"


###############################################################################
#------------------rules
###############################################################################

rule all:
  input:
    # Figure 1
    "results/current/Figures/Figure_chipseq.pdf",
    # Figure 2
    "results/current/Figures/Figure_peakgeneannotation.pdf",
    "results/current/integration/deeptools/GR_chipseq_UPandDOWN_libnorm_30kb_1000b_1000a_heatmap.svg",
    # Figure 3
    "results/current/Figures/Figure_abcresults.pdf",
    # Figure 4
    "results/current/Figures/Figure_GLMs.pdf",
    # Figure 5 (A+B)
    "results/current/Figures/Figure_stats.pdf",
    # Suppl. Figure 1
    "results/current/Figures/Figure_rnaseq.pdf",
    # Suppl. Figure 2
    "results/current/Figures/Figure_motifcorrelations.pdf",
    # Suppl. Figure 4
    "results/current/Figures/Figure_supplemental_GLMs.pdf",
    # Suppl. Figure 5
    "results/current/Figures/Suppl_Figure_statmotifs.pdf",

    # tracks for visualizations with IGV
    "results/current/ChIP/GR_2020/bw/DexLPS_chipseq_GR_PE_merged_GRCm38_cpm.bw",
    "results/current/ChIP/H3K27ac/bw/DexLPS_histone_H3K27ac_PE_merged_GRCm38_libnorm.bw",
    "results/current/ChIP/H3K27ac/bw/LPS_histone_H3K27ac_PE_merged_GRCm38_libnorm.bw",
    "results/current/atacseq/bw/merged_LPS_GRCm38_libnorm.bw",
    "results/current/atacseq/bw/merged_DexLPS_GRCm38_libnorm.bw",


rule download_JASPAR:
    input:
    output:
        "data/current/motifs/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt"
    params:
        outdir="data/current/motifs/",
        url="https://jaspar.genereg.net/download/data/2022/CORE/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt"
    shell:
        "wget {params.url} -P {params.outdir}"

rule make_bowtie_index:
    input:
        fa=config['genome']
    output:
        expand("data/current/bowtie_index/{basename}.1.ebwt", basename=config['bowtie']['basename'])
        # the resulting indeces are faster versions of .fa files, specific to bowtie1
    params:
        basename=config['bowtie']['basename'],
        outdir="data/current/bowtie_index",
    conda:
        "src/conda/samtools_19.yml"
    shell:
        "bowtie-build -f {input[fa]} {params.outdir}/{params.basename}"

rule get_reference_files:
    input:
    output:
        "data/current/mm10.blacklist.bed.gz",
        "data/current/Mus_musculus.GRCm38.100.gtf.gz",
        "data/current/mm9ToMm10.over.chain",
        "data/current/gencode_annotations/gencode.vM1.annotation.gtf.gz",
        "data/current/gencode_annotations/gencode.vM23.annotation.gene.gtf"
    params:
        outdir="data/current"
    shell:
        """
        wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm10-mouse/mm10.blacklist.bed.gz -P {params.outdir} && \
        wget http://ftp.ensembl.org/pub/release-100/gtf/mus_musculus/Mus_musculus.GRCm38.100.gtf.gz -P {params.outdir} && \
        wget http://hgdownload.cse.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz -P {params.outdir} && \
        gunzip data/current/mm9ToMm10.over.chain.gz && \

        mkdir -p "data/current/gencode_annotations" && \
        wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M1/gencode.vM1.annotation.gtf.gz -P {params.outdir}/gencode_annotations/ && \
        wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.annotation.gtf.gz -P {params.outdir}/gencode_annotations/ && \
        zcat data/current/gencode_annotations/gencode.vM23.annotation.gtf.gz | grep -P '\tgene\t' > data/current/gencode_annotations/gencode.vM23.annotation.gene.gtf
        """  

rule get_remap_data:
    input:
    output:
        fullbed = "data/current/remap/remap2022_all_macs2_mm10_v1_0.bed",
        separatedfulltsv = "data/current/remap/remap2022_all_macs2_mm10_v1_0_separated.tsv",
        filteredtsv = "data/current/remap/remap2022_all_macs2_mm10_v1_0_separated_ftr.tsv",
        filteredbed = "data/current/remap/remap2022_all_macs2_mm10_v1_0_separated_ftr.bed"
    params:
        downloadlink="https://remap.univ-amu.fr/storage/remap2022/mm10/MACS2/remap2022_all_macs2_mm10_v1_0.bed.gz",
        outdir="data/current/remap" 
    shell:
        """
        mkdir -p {params.outdir} && \
        wget {params.downloadlink} -P {params.outdir} && \
        gunzip -c {params.outdir}/remap2022_all_macs2_mm10_v1_0.bed.gz > {output.fullbed} && \
        sed -E "s/(\\w+)\\.(\\w+)\\.(\\w+)/\\1\\t\\2\\t\\3/g" {output.fullbed} > {output.separatedfulltsv} && \
        awk -F'\t' '$6~/BMDM|macrophage/' {output.separatedfulltsv} > {output.filteredtsv} && \
        awk '{{print $1"\t"$2"\t"$3"\t"$4","$5","$6}}' FS='\t' {output.filteredtsv} > {output.filteredbed}
        """

rule gunzip_blacklist:
    input:
        "data/current/mm10.blacklist.bed.gz"
    output:
        "data/current/mm10.blacklist.bed"
    shell:
        "gunzip -c {input} > {output}" # -c keeps the original file unchanged

rule scan_nr3c1_genomewide:
    input:
        motif="data/current/motifs/custom/nr3c1_simplified_{sitelength}.motif",
        fa= "data/current/genome/GRCm38.primary_assembly.standchr.fa"
    output:
        "results/current/homer/nr3c1.{sitelength}matches.mm10.bed"
    shell:
        """
        mkdir -p $(dirname {output}) && \
        scanMotifGenomeWide.pl {input.motif} {input.fa} -bed > {output}
        """
  