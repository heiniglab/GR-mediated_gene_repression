include: "Snakefile_utils.smk"

rule merge_GR2020:
    input: 
        "data/current/chipseq/GR_2020/raw_data/{sampledir}/{prefix}_L001_R{suffix}.fastq.gz",
        "data/current/chipseq/GR_2020/raw_data/{sampledir}/{prefix}_L002_R{suffix}.fastq.gz"
    output:
        "data/current/chipseq/GR_2020/fastq_merged/{sampledir}/{prefix}_R{suffix}.fastq.gz"
    shell:
        """
        mkdir -p $(dirname {output}) && \
        cat {input} > {output}
        """

rule symlink_GR2020:
    input:
        # Input controls
        "data/current/chipseq/GR_2020/input_ctrls/GAR0517/fastq/GAR0517_BC7FEMANXX_AGTCAA_L007_R1_001.fastq.gz",
        "data/current/chipseq/GR_2020/input_ctrls/GAR0517/fastq/GAR0517_BC7FEMANXX_AGTCAA_L007_R2_001.fastq.gz",
        "data/current/chipseq/GR_2020/input_ctrls/Sample_MUC9117/MUC9117_R1_merged.fastq.gz",
        "data/current/chipseq/GR_2020/input_ctrls/Sample_MUC9117/MUC9117_R2_merged.fastq.gz",
        "data/current/chipseq/GR_2020/input_ctrls/Sample_MUC9118/MUC9118_R1_merged.fastq.gz",
        "data/current/chipseq/GR_2020/input_ctrls/Sample_MUC9118/MUC9118_R2_merged.fastq.gz",
        "data/current/chipseq/GR_2020/input_ctrls/Sample_GAR1531/GAR1531_S13_L002_R1_001.fastq.gz",
        "data/current/chipseq/GR_2020/input_ctrls/Sample_GAR1531/GAR1531_S13_L002_R2_001.fastq.gz",
        # Samples
        "data/current/chipseq/GR_2020/fastq_merged/Sample_MUC20387/MUC20387_S6_R1_001.fastq.gz",
        "data/current/chipseq/GR_2020/fastq_merged/Sample_MUC20387/MUC20387_S6_R2_001.fastq.gz",
        "data/current/chipseq/GR_2020/fastq_merged/Sample_MUC20388/MUC20388_S7_R1_001.fastq.gz",
        "data/current/chipseq/GR_2020/fastq_merged/Sample_MUC20388/MUC20388_S7_R2_001.fastq.gz"
    output:
        expand("results/current/ChIP/GR_2020/fastq/raw/DexLPS_chipseq_GR_rep{rep}_R{dir}.fastq.gz", rep=[1,2], dir=[1,2] ),
        expand("results/current/ChIP/GR_2020/fastq/raw/input_rep{rep}_R{dir}.fastq.gz", rep=[1,2,3,4], dir=[1,2] )  
    params:
        indir_ctrl="data/current/chipseq/GR_2020/input_ctrls",
        indir_samples="data/current/chipseq/GR_2020/fastq_merged",
        outdir="results/current/ChIP/GR_2020/fastq/raw"
    shell:
        """
        mkdir -p {params.outdir} && \
        sh src/scripts/symlink_GR_2020.sh {params.indir_ctrl} {params.indir_samples} {params.outdir}
        """

rule symlink_histone:
    input:
        "data/current/chipseq/H3K27ac/fastq/GAR0814_S7_L002_R1_001.fastq.gz",
        "data/current/chipseq/H3K27ac/fastq/GAR0814_S7_L002_R2_001.fastq.gz",
        "data/current/chipseq/H3K27ac/fastq/GAR0815_S8_L002_R1_001.fastq.gz",
        "data/current/chipseq/H3K27ac/fastq/GAR0815_S8_L002_R2_001.fastq.gz",
        "data/current/chipseq/H3K27ac/fastq/GAR0816_S9_L002_R1_001.fastq.gz",
        "data/current/chipseq/H3K27ac/fastq/GAR0816_S9_L002_R2_001.fastq.gz",
        "data/current/chipseq/H3K27ac/fastq/GAR0823_S10_L002_R1_001.fastq.gz",
        "data/current/chipseq/H3K27ac/fastq/GAR0823_S10_L002_R2_001.fastq.gz"
    output:
        "results/current/ChIP/H3K27ac/fastq/raw/LPS_histone_H3K27ac_rep1_R1.fastq.gz",
        "results/current/ChIP/H3K27ac/fastq/raw/LPS_histone_H3K27ac_rep1_R2.fastq.gz",
        "results/current/ChIP/H3K27ac/fastq/raw/LPS_histone_H3K27ac_rep2_R1.fastq.gz",
        "results/current/ChIP/H3K27ac/fastq/raw/LPS_histone_H3K27ac_rep2_R2.fastq.gz",
        "results/current/ChIP/H3K27ac/fastq/raw/DexLPS_histone_H3K27ac_rep1_R1.fastq.gz",
        "results/current/ChIP/H3K27ac/fastq/raw/DexLPS_histone_H3K27ac_rep1_R2.fastq.gz",
        "results/current/ChIP/H3K27ac/fastq/raw/DexLPS_histone_H3K27ac_rep2_R1.fastq.gz",
        "results/current/ChIP/H3K27ac/fastq/raw/DexLPS_histone_H3K27ac_rep2_R2.fastq.gz"
    params:
        indir="data/current/chipseq/H3K27ac/fastq",
        outdir="results/current/ChIP/H3K27ac/fastq/raw"
    shell:
        """
        mkdir -p {params.outdir} && \
        sh src/scripts/symlink_histone_samples.sh {params.indir} {params.outdir}
        """

rule symlink_test:
    input:
        "data/current/chipseq/H3K27ac/fastq/GAR0814_S7_L002_R1_001.fastq.gz",
        "data/current/chipseq/H3K27ac/fastq/GAR0814_S7_L002_R2_001.fastq.gz",
        "data/current/chipseq/H3K27ac/fastq/GAR0815_S8_L002_R1_001.fastq.gz",
        "data/current/chipseq/H3K27ac/fastq/GAR0815_S8_L002_R2_001.fastq.gz",
        "data/current/chipseq/H3K27ac/fastq/GAR0816_S9_L002_R1_001.fastq.gz",
        "data/current/chipseq/H3K27ac/fastq/GAR0816_S9_L002_R2_001.fastq.gz",
        "data/current/chipseq/H3K27ac/fastq/GAR0823_S10_L002_R1_001.fastq.gz",
        "data/current/chipseq/H3K27ac/fastq/GAR0823_S10_L002_R2_001.fastq.gz"
    output:
        "results/current/ChIP/H3K27ac/fastq/raw/test"
    params:
        indir="data/current/chipseq/H3K27ac/fastq",
        outdir="results/current/ChIP/H3K27ac/fastq/raw"
    shell:
        """
        mkdir -p {params.outdir} && \
        sh src/scripts/symlink_test.sh {params.indir} {params.outdir}
        """

#--------------------------------------------------------------
#---processing
#--------------------------------------------------------------
rule trim_adapter_sequences:
    input:
        fastq="results/current/ChIP/{exp}/fastq/raw/{treatment}_chipseq_{target}.fastq.gz"
    output:
        "results/current/ChIP/{exp}/fastq/trimmed/{treatment}_chipseq_{target}_trimmed.fastq.gz",
    params:
        minimum_overlap = config['fastq_trim_adapters']['minimum_overlap'], #minimum overlap allowed
        minimum_fragment_length = config['fastq_trim_adapters']['minimum_fragment_length'], #keep fragments longer than this length
        maximum_error_rate = config['fastq_trim_adapters']['maximum_error_rate'], #maximum error rate allowed for adapter searches
    message:
        "Trimming adapters from fastq reads..."
    shell:
        "mkdir -p results/current/ChIP/{wildcards.exp}/fastq/trimmed; \
        cutadapt -m {params.minimum_fragment_length} \
        -a AGATCGGAAGAGCACACGTCT \
        -o {output} {input.fastq}"

#--------------------------------------------------------------
#------------multiQC
#--------------------------------------------------------------

# run multiqc after the mapping and fastqc of all samples has finished, so that gets included in the output   
rule multiqc:
    input:
        expand("results/current/ChIP/H3K27ac/bam/{condition}_histone_{factor}_rep{rep}_{basename}_algnd.bam",  condition=['LPS','DexLPS'], factor=['H3K27ac'], rep=[1,2], basename=config['bowtie']['basename']),
        expand("results/current/ChIP/H3K27ac/QC/{condition}_histone_{factor}_R{r}_fastqc.html", condition=['LPS','DexLPS'], factor=['H3K27ac'], r=[1,2])
    output:
        "results/current/ChIP/histone/multiqc_report.html"
    message:
        "Running MultiQC on histone data"
    shell:
        "multiqc results/current/logs results/current/ChIP/H3K27ac/QC --outdir results/current/ChIP/H3K27ac --force"  

#------------------------------------------
#-- alignment, BAM filtering and sorting
#------------------------------------------

rule fastq_align:
    input:
        fastq="results/current/ChIP/{exp}/fastq/trimmed/{treatment}_chipseq_{target}_trimmed.fastq.gz",
        ebwt="data/current/bowtie_index/{basename}.1.ebwt", # check if the index files have been created
    output:
        "results/current/ChIP/{exp}/bam/{treatment}_chipseq_{target}_{basename}_algnd.bam",
    params:
        index=expand("data/current/bowtie_index/{basename}", basename=config['bowtie']['basename']),
        nworkers = config['fastq_align']['cpus'], #numer of parallel processes
        alignments_reported = 1, #how many valid alignments reported
        multialignments_allowed = 1, #suppress all alignments if more than these reportable alignments exist
        mismatches_allowed = config['fastq_align']['mismatches_allowed'], #mismatches allowed in alignments
    conda:
        "../conda/samtools_19.yml"
    shell:
        "mkdir -p results/current/logs && \
        mkdir -p results/current/ChIP/{wildcards.exp}/bam && \
        gzip -cd {input[fastq]} | bowtie -S -p {params.nworkers} --chunkmbs 512 \
        -k {params.alignments_reported} -m {params.multialignments_allowed} \
        -v {params.mismatches_allowed} {params.index} - | \
        samtools view -F 4 -Sbo {output} -"
 
# Align ChIP-input samples using bowtie1
rule PE_fastq_align:
    input:
        fa1="{path1}/fastq/raw/{sample}_R1.fastq.gz",
        fa2="{path1}/fastq/raw/{sample}_R2.fastq.gz",
        bowtie_index=expand("data/current/bowtie_index/{basename}.1.ebwt", basename=config['bowtie']['basename'])
    output:
        "{path1}/bam/{sample}_PE_{basename}_algnd.bam",
        # "results/current/ChIP/H3K27ac/bam/LPS_histone_H3K27ac_rep1_PE_GRCm38_algnd.bam"
    params:
        index= expand("data/current/bowtie_index/{basename}",basename=config['bowtie']['basename']),
        nworkers = config['fastq_align']['cpus'], #numer of parallel processes
        alignments_reported = 1, #how many valid alignments reported
        multialignments_allowed = 1, #suppress all alignments if more than these reportable alignments exist
        mismatches_allowed = config['fastq_align']['mismatches_allowed'], #mismatches allowed in alignments
    message:
        "Aligning PE sample...",
    conda:
        "../conda/samtools_19.yml"
    shell:
        "mkdir -p results/current/logs && \
        mkdir -p $(dirname {output}) && \
        bowtie -S -p {params.nworkers} --chunkmbs 512 \
        -k {params.alignments_reported} -m {params.multialignments_allowed} \
        -v {params.mismatches_allowed} {params.index} \
        -1 {input[fa1]} -2 {input[fa2]} | samtools view -F 4 -Sbo {output} -"

# rules for sorting, blacklist filtering and deduplication of the resulting bam files are defined within `Snakefile_utils.smk`

rule run_macs2_narrow:
    input:
        sample="results/current/ChIP/{exp}/bam/{treatment}_chipseq_{target}_sorted_ftr_dedup.bam",
        control="results/current/ChIP/{exp}/bam/input_sorted_ftr_dedup.bam"
    output:
        "results/current/ChIP/{exp}/macs2/{treatment}_chipseq_{target}_q0.05_peaks.narrowPeak",
        "results/current/ChIP/{exp}/macs2/{treatment}_chipseq_{target}_q0.05_summits.bed"
        #"results/current/ChIP/cofactors2013/macs2/DexLPS_chipseq_cJun_GRCm38_q0.05_summits.bed"
    params:
        name = "{treatment}_chipseq_{target}_q0.05",
        tempdir = "/localscratch/barbara.hoellbacher/tmp",
        outdir = "results/current/ChIP/{exp}/macs2"
    shell:
        """
        mkdir -p {params.tempdir} && \
        mkdir -p {params.outdir} && \
        macs2 callpeak -t {input.sample} -c {input.control}\
        --tempdir {params.tempdir} --qvalue 0.05 \
        --keep-dup all --bdg --SPMR \
        -g mm -f BAM --outdir {params.outdir} -n {params.name}
        """

rule run_macs2_broad:
    input:
        sample="results/current/ChIP/{target}/bam/{treatment}_histone_{target}_PE_merged_GRCm38_sorted_ftr_dedup.bam"
    output:
        "results/current/ChIP/{target}/macs2/{treatment}_histone_{target}_q0.05_PE_merged_peaks.broadPeak"
    params:
        name = "{treatment}_histone_{target}_q0.05_PE_merged",
        tempdir = "/localscratch/barbara.hoellbacher/tmp",
        outdir = "results/current/ChIP/{target}/macs2"
    shell:
        """
        mkdir -p {params.tempdir} && \
        mkdir -p {params.outdir} && \
        macs2 callpeak -t {input.sample} --broad \
        --tempdir {params.tempdir} --qvalue 0.05 \
        --keep-dup all --bdg --SPMR \
        -g mm -f BAM --outdir {params.outdir} -n {params.name}
        """

# Run MACS2 in PE mode with relaxed threshold for IDR
rule run_macs2_GR2020:
    input:
        "results/current/ChIP/GR_2020/bam/DexLPS_chipseq_GR_rep{rep}_PE_{basename}_sorted_ftr_dedup.bam",
        "results/current/ChIP/GR_2020/bam/input_rep1_PE_{basename}_sorted_ftr_dedup.bam",
        "results/current/ChIP/GR_2020/bam/input_rep2_PE_{basename}_sorted_ftr_dedup.bam",
        "results/current/ChIP/GR_2020/bam/input_rep3_PE_{basename}_sorted_ftr_dedup.bam",
        "results/current/ChIP/GR_2020/bam/input_rep4_PE_{basename}_sorted_ftr_dedup.bam",
    output:
        "results/current/ChIP/GR_2020/macs2/DexLPS_chipseq_GR_rep{rep}_PE_{basename}_p0.1_peaks.narrowPeak",
        "results/current/ChIP/GR_2020/macs2/DexLPS_chipseq_GR_rep{rep}_PE_{basename}_p0.1_summits.bed"
    params:
        name = "DexLPS_chipseq_GR_rep{rep}_PE_{basename}_p0.1",
        outdir = "results/current/ChIP/GR_2020/macs2",
        tempdir = "/localscratch/barbara.hoellbacher/tmp"
    message: "Running MACS2 with relaxed threshold for IDR..."
    shell:
        "mkdir -p $(dirname {output[0]}) && \
        mkdir -p {params.tempdir} && \
        macs2 callpeak -t {input[0]} -c {input[1]} {input[2]} {input[3]} {input[4]}\
        --tempdir {params.tempdir} --pvalue 0.1 \
        --keep-dup 1 --bdg -g mm -f BAMPE --outdir {params.outdir} -n {params.name}"

# run IDR on replicates
rule GRChIP_run_idr_for_replicates:
    input:
        "results/current/ChIP/GR_2020/macs2/DexLPS_chipseq_GR_rep1_PE_{basename}_p0.1_peaks.narrowPeak",
        "results/current/ChIP/GR_2020/macs2/DexLPS_chipseq_GR_rep2_PE_{basename}_p0.1_peaks.narrowPeak",
    output:
        "results/current/ChIP/GR_2020/idr/DexLPS_chipseq_GR_{basename}_p0.1_idr0-05"
    message: "Running IDR on replicates"
    shell:
        "mkdir -p $(dirname {output}) && \
        idr --samples {input} --input-file-type narrowPeak --rank signal.value --idr-threshold 0.05 --output-file {output} --plot"

# merge BAMs to then create ONE bigwig file as input to the NN
rule merge_BAMs_from_replicates:
    input:
        "results/current/ChIP/{exp}/bam/{sample}_rep1_PE_{basename}_sorted_ftr_dedup.bam",
        "results/current/ChIP/{exp}/bam/{sample}_rep2_PE_{basename}_sorted_ftr_dedup.bam",
    output:
        "results/current/ChIP/{exp}/bam/{sample}_PE_merged_{basename}_sorted_ftr_dedup.bam",
    conda:
        "../conda/samtools_19.yml"
    shell:
        "samtools merge {output} {input[0]} {input[1]}"

rule get_adjusted_libsize_atac:
    input:
        peaks1="results/current/atacseq/macs2/DexLPS_peaks.narrowPeak.sorted",
        peaks2="results/current/atacseq/macs2/LPS_peaks.narrowPeak.sorted",
        bam1="results/current/atacseq/bam/merged_DexLPS_GRCm38_sorted_ftr_dedup.bam",
        bam2="results/current/atacseq/bam/merged_LPS_GRCm38_sorted_ftr_dedup.bam"
    output:
        peakuniverse="results/current/atacseq/macs2/atac_universe.narrowPeak.sorted",
        counts="results/current/atacseq/counts",
    shell:
        """
        cat {input.peaks1} {input.peaks2} | bedtools sort -i stdin | bedtools merge -i stdin > {output.peakuniverse} && \
        echo {input.bam1} > {output.counts} && \
        samtools view -c -L {output.peakuniverse} {input.bam1} >> {output.counts} && \
        echo {input.bam2} >> {output.counts} && \
        samtools view -c -L {output.peakuniverse} {input.bam2} >> {output.counts}
        """
  
rule get_adjusted_libsize_h3k27ac:
    input:
        peaks1="results/current/ChIP/H3K27ac/macs2/DexLPS_histone_H3K27ac_q0.05_PE_merged_peaks.broadPeak",
        peaks2="results/current/ChIP/H3K27ac/macs2/LPS_histone_H3K27ac_q0.05_PE_merged_peaks.broadPeak",
        bam1="results/current/ChIP/H3K27ac/bam/DexLPS_histone_H3K27ac_PE_merged_GRCm38_sorted_ftr_dedup.bam",
        bam2="results/current/ChIP/H3K27ac/bam/LPS_histone_H3K27ac_PE_merged_GRCm38_sorted_ftr_dedup.bam"
    output:
        peakuniverse="results/current/ChIP/H3K27ac/macs2/h3k27ac_universe.broadPeak",
        counts="results/current/ChIP/H3K27ac/counts"
    shell:
        """
        cat {input.peaks1} {input.peaks2} | bedtools sort -i stdin | bedtools merge -i stdin > {output.peakuniverse} && \
        echo {input.bam1} > {output.counts} && \
        samtools view -c -L {output.peakuniverse} {input.bam1} >> {output.counts} && \
        echo {input.bam2} >> {output.counts} && \
        samtools view -c -L {output.peakuniverse} {input.bam2} >> {output.counts}
        """
