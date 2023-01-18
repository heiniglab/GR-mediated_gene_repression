echo "Creating Profileplot"

set -eo pipefail

input=$1
bw_pos=$2
bw_neg=$3
output=$4
outdir=$(dirname $4)

awk 'NR > 1 {{print $3,$4,$5,$6}}' FS='\t' OFS='\t' ${input} \
> ${outdir}/fimo_motifmatched.bed

computeMatrix reference-point --referencePoint center \
-S ${bw_pos} ${bw_neg} \
-R ${outdir}/fimo_motifmatched.bed \
-a 100 -b 100 -o  ${outdir}/computeMatrix_mergedfilteredbws

plotProfile --matrixFile  ${outdir}/computeMatrix_mergedfilteredbws  \
--samplesLabel pos_strand neg_strand \
--plotTitle "5' end distro w.r.t. NR3C1 motif" \
--refPointLabel motifcenter \
--perGroup --colors blue olive -o ${output}
