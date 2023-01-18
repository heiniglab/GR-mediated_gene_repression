# how these trimmed files are generated depends on the experiment
rule generate_fastqc_report:
    input:
        "{path1}/fastq/trimmed/{path2}_trimmed.fastq.gz"
    output:
        "{path1}/QC/{path2}_trimmed_fastqc.html"
    message:
        "Generate fastqc reports"
    shell:
        "mkdir -p {wildcards.path1}/QC; \
        fastqc -o {wildcards.path1}/QC -f fastq {input}"

rule sort_raw_BAM_and_blacklist_filter:
    input:
        bam="{file}_algnd.bam",
        blacklist="data/current/mm10.blacklist.bed.gz"
    output:
        "{file}_sorted_ftr.bam"
    conda:
        "../conda/samtools_19.yml"
    threads: 4
    message:
        "Sorting the raw bam files and removing blacklisted regions"
    shell:
        "mkdir -p results/current/tmp && \
        samtools sort {input.bam} -T results/current/tmp/ --threads {threads} | \
        bedtools intersect -v -abam stdin -b {input.blacklist} > {output}"

rule filter_duplicates:
    input:
        "{path1}/bam/{path2}_sorted_ftr.bam"
    output:
        bam = "{path1}/bam/{path2}_sorted_ftr_dedup.bam",
        metric = "{path1}/stats/{path2}.picardDupMetrics.txt"
    shell:
        "mkdir -p {wildcards.path1}/stats && \
        picard MarkDuplicates INPUT={input} OUTPUT={output.bam} REMOVE_DUPLICATES=true METRICS_FILE={output.metric} VALIDATION_STRINGENCY=LENIENT PROGRAM_RECORD_ID='null'"

rule make_bw_from_bam:
    input:
        bam="{path1}/bam/{path2}_GRCm38_sorted_ftr_dedup.bam",
        bai="{path1}/bam/{path2}_GRCm38_sorted_ftr_dedup.bai"
    output:
        "{path1}/bw/{path2}_GRCm38.bw"
    params:
        processors=config['bamCoverage']['cpus']
    shell:
        "bamCoverage -b {input.bam} --normalizeUsing None \
        --binSize 1 -p {params.processors} --outFileName {output[0]}" 

rule make_forwardbw_from_bam:
    input:
        bam="{path1}/bam/{path2}_GRCm38_sorted_ftr_dedup.bam",
        bai="{path1}/bam/{path2}_GRCm38_sorted_ftr_dedup.bai"
    output:
        "{path1}/bw/{path2}_GRCm38_forward.bw"
    params:
        processors=config['bamCoverage']['cpus']
    shell:
        "bamCoverage -b {input.bam} --normalizeUsing None --filterRNAstrand forward \
        --binSize 1 -p {params.processors} --outFileName {output[0]}" 
        
rule make_reversebw_from_bam:
    input:
        bam="{path1}/bam/{path2}_GRCm38_sorted_ftr_dedup.bam",
        bai="{path1}/bam/{path2}_GRCm38_sorted_ftr_dedup.bai"
    output:
        "{path1}/bw/{path2}_GRCm38_reverse.bw"
    params:
        processors=config['bamCoverage']['cpus']
    shell:
        "bamCoverage -b {input.bam} --normalizeUsing None --filterRNAstrand reverse \
        --binSize 1 -p {params.processors} --outFileName {output[0]}" 

rule make_RPGC_bw_from_bam:
    input:
        bam="{path1}/bam/{path2}_GRCm38_sorted_ftr_dedup.bam",
        bai="{path1}/bam/{path2}_GRCm38_sorted_ftr_dedup.bai"
    output:
        "{path1}/bw/{path2}_GRCm38_rpgc.bw"
    params:
        processors=config['bamCoverage']['cpus']
    shell:
    # the specified effective genomesize is a good estimate for read length of 100bp, see https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html
        "bamCoverage -b {input.bam} --normalizeUsing RPGC --effectiveGenomeSize 2467481108 \
        --binSize 1 -p {params.processors} --outFileName {output[0]}" # normalization: RPGC = reads per genomic content (1x normalization) scales for 1x average coverage

rule make_CPM_bw_from_bam:
    input:
        bam="{path1}/bam/{path2}_GRCm38_sorted_ftr_dedup.bam",
        bai="{path1}/bam/{path2}_GRCm38_sorted_ftr_dedup.bai"
    output:
        "{path1}/bw/{path2}_GRCm38_cpm.bw"
    params:
        processors=config['bamCoverage']['cpus']
    shell:
        "bamCoverage -b {input.bam} --normalizeUsing CPM \
        --binSize 1 -p {params.processors} --outFileName {output[0]}" # normalization: CPM: Counts Per Million mapped reads

# generate bed files from IDR peaks for motif analysis
rule generate_bed_from_idr:
    input:
        "{path1}/idr/{path2}_idr0-05"
    output:
        "{path1}/idr/{path2}_idr0-05.bed"
    message: "Parse genomic locations from narrowPeak files"
    shell:
        "cut -f 1,2,3,4,5 {input} > {output}"

# sort narrowPeak files by -log10(p-value)
rule sort_narrowPeak_bypval:
    input:
        "{file}_peaks.narrowPeak"
    output:
        "{file}_peaks_sorted.narrowPeak"
    message: "Sorting narrowPeak files by p-value for IDR"
    shell:
        "sort -k8,8nr {input} > {output}"

rule read_summit_from_IDR_peaks:
    input:
        "{file}_idr0-05"
    output:
        "{file}_idr0-05_summit"
    message: "Write summit location of peaks passing threshold"
    shell:
        "awk '{{print $1,$2+$10,$2+$10+1,$4,$5}}' FS='\\t' OFS='\t' {input} > {output}"

rule get_chromsize_from_fa:
    input:
        fa=config['genome'],
    output:
        chromsizes=expand("results/current/genomesize/{basename}.chromsizes", basename=config['bowtie']['basename'])
    conda:
        "../conda/samtools_19.yml"
    params:
    shell:
        "mkdir -p results/current/genomesize/; \
        samtools faidx {input[fa]} && \
        cut -f1,2 {input[fa]}.fai > {output[chromsizes]}"

rule get_chromsize_bed:
    input:
        expand("results/current/genomesize/{basename}.chromsizes", basename=config['bowtie']['basename'])
    output:
         expand("results/current/genomesize/{basename}.chromsizes.bed", basename=config['bowtie']['basename'])
    shell:
        "awk '{{print $1,0,$2}}' FS='\\t' OFS='\t' {input} > {output}"

rule create_bed_interval_from_summit:
    input:
        summit="{file}_idr0-05_summit",
        chromsizes=expand("results/current/genomesize/{basename}.chromsizes", basename=config['bowtie']['basename'])
    output:
        bed="{file}_idr0-05_summitinterval_slop{slop}_unmerged.bed",
    message:
        "Create bed interval from summit"
    shell:
        """
        bedtools slop -i {input[summit]} -g {input[chromsizes]} -b {wildcards.slop} | \
        sort -k1,1 -k2,2n - > {output[bed]}
        """

rule merge_overlapping_bed_interval:
    input:
        bed="{file}_idr0-05_summitinterval_slop{slop}_unmerged.bed"
    output:
        bed="{file}_idr0-05_summitinterval_slop{slop}_merged.bed"
    message:
        "Merge overlapping bed intervals"
    shell:
        """
        bedtools merge -i {input.bed} -c 4 -o collapse -delim "|" > {output.bed}
        """

rule bam_index:
    input:
        "{file}.bam"
    output:
        "{file}.bai"
    conda:
        "../conda/samtools_19.yml"
    shell:
        "samtools index {input} && mv {input}.bai {output}"
