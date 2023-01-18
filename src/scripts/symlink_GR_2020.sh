#!/bin/bash
set -eo pipefail

# for the input samples
indir_ctrl=$1
indir_samples=$2
outdir=$3

ln -s "$(realpath $indir_ctrl/GAR0517/fastq/GAR0517_BC7FEMANXX_AGTCAA_L007_R1_001.fastq.gz)" $outdir/input_rep1_R1.fastq.gz
ln -s "$(realpath $indir_ctrl/GAR0517/fastq/GAR0517_BC7FEMANXX_AGTCAA_L007_R2_001.fastq.gz)" $outdir/input_rep1_R2.fastq.gz
ln -s "$(realpath $indir_ctrl/Sample_MUC9117/MUC9117_R1_merged.fastq.gz)" $outdir/input_rep2_R1.fastq.gz
ln -s "$(realpath $indir_ctrl/Sample_MUC9117/MUC9117_R2_merged.fastq.gz)" $outdir/input_rep2_R2.fastq.gz
ln -s "$(realpath $indir_ctrl/Sample_MUC9118/MUC9118_R1_merged.fastq.gz)" $outdir/input_rep3_R1.fastq.gz
ln -s "$(realpath $indir_ctrl/Sample_MUC9118/MUC9118_R2_merged.fastq.gz)" $outdir/input_rep3_R2.fastq.gz
ln -s "$(realpath $indir_ctrl/Sample_GAR1531/GAR1531_S13_L002_R1_001.fastq.gz)" $outdir/input_rep4_R1.fastq.gz
ln -s "$(realpath $indir_ctrl/Sample_GAR1531/GAR1531_S13_L002_R2_001.fastq.gz)" $outdir/input_rep4_R2.fastq.gz

# for the GR 2020 GR samples
ln -s "$(realpath $indir_samples/Sample_MUC20387/MUC20387_S6_R1_001.fastq.gz)" $outdir/DexLPS_chipseq_GR_rep1_R1.fastq.gz
ln -s "$(realpath $indir_samples/Sample_MUC20387/MUC20387_S6_R2_001.fastq.gz)" $outdir/DexLPS_chipseq_GR_rep1_R2.fastq.gz
ln -s "$(realpath $indir_samples/Sample_MUC20388/MUC20388_S7_R1_001.fastq.gz)" $outdir/DexLPS_chipseq_GR_rep2_R1.fastq.gz
ln -s "$(realpath $indir_samples/Sample_MUC20388/MUC20388_S7_R2_001.fastq.gz)" $outdir/DexLPS_chipseq_GR_rep2_R2.fastq.gz