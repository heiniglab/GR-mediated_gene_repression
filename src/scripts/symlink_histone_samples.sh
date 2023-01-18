#!/bin/bash
set -eo pipefail

# for the input samples
indir=$1
outdir=$2

ln -s "$(realpath $indir/GAR0814_S7_L002_R1_001.fastq.gz)" $outdir/LPS_histone_H3K27ac_rep1_R1.fastq.gz
ln -s "$(realpath $indir/GAR0814_S7_L002_R2_001.fastq.gz)" $outdir/LPS_histone_H3K27ac_rep1_R2.fastq.gz
ln -s "$(realpath $indir/GAR0815_S8_L002_R1_001.fastq.gz)" $outdir/LPS_histone_H3K27ac_rep2_R1.fastq.gz
ln -s "$(realpath $indir/GAR0815_S8_L002_R2_001.fastq.gz)" $outdir/LPS_histone_H3K27ac_rep2_R2.fastq.gz
ln -s "$(realpath $indir/GAR0816_S9_L002_R1_001.fastq.gz)" $outdir/DexLPS_histone_H3K27ac_rep1_R1.fastq.gz
ln -s "$(realpath $indir/GAR0816_S9_L002_R2_001.fastq.gz)" $outdir/DexLPS_histone_H3K27ac_rep1_R2.fastq.gz
ln -s "$(realpath $indir/GAR0823_S10_L002_R1_001.fastq.gz)" $outdir/DexLPS_histone_H3K27ac_rep2_R1.fastq.gz
ln -s "$(realpath $indir/GAR0823_S10_L002_R2_001.fastq.gz)" $outdir/DexLPS_histone_H3K27ac_rep2_R2.fastq.gz

