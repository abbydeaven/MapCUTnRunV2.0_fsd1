#!/bin/bash
#SBATCH --partition=batch
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ad45368@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100gb
#SBATCH --time=8:00:00
#SBATCH --output=../MappingOutput/logs/%x.out
#SBATCH --error=../MappingOutput/logs/%x.err


cd $SLURM_SUBMIT_DIR
##Directory information+variables:
dir=/lustre2/scratch/ad45368/ChIPseq/FinalH3k27me3_ChIP
outdir="/lustre2/scratch/ad45368/ChIPseq/FinalH3k27me3_ChIP"
bamdir="${outdir}/bamFiles"
bwdir="${outdir}/bigWig"
PeakDir="${outdir}/Peaks"
Motfis="${outdir}/Motifs"

## Input control and rep information, using name:
dir="/lustre2/scratch/ad45368/ChIPseq/PeritheciaChIPAnalysis/2025_ChIP_Reanalysis/MappingOutput/bamFiles"
#Control1="143_144_ChIP_WT_gfp_trap_Rep1"
dpf3_input_1='145-43_ChIP_WT_3dpf_input_Rep1'
dpf3_input_2='153_43_ChIP_WT_dpf3_Input_Rep1_S43_L002'

dpf3_H3K27me3_1='145_44_ChIP_WT_3dpf_H3K27me3_Rep1'
dpf3_H3K27me3_2='146_111_ChIP_WT_3dpf_H3K27me3_Rep2'
dpf3_H3K27me3_3='152_35_ChIP_WT_3dpf_P_H3K27me3_Rep6'

dpf6_H3K27me3_1='150-28_ChIP_WT_P_6pf_H3K27me3_Rep2'
dpf6_H3K27me3_2='151_128_ChIP_WT_P_6dpf_H3K27me3_Rep3'
dpf6_H3K27me3_3='152_38_ChIP_WT_6dpf_P_H3K27me3_Rep5'

Mycelia_H3K27me3_1='150-31_ChIP_WT_M_H3K27me3_Rep1'
Mycelia_H3K27me3_2='151-130_ChIP_WT_M_H3K27me3_Rep2'
## Use bedtools to make a consensus peakset of peaks with 80% overlap in all samples
ml BEDTools/2.31.1-GCC-13.3.0

bedtools intersect -a $PeakDir/${dpf3_input_1}_peaks.broadPeak -b ${PeakDir}/${dpf3_input_2}_peaks.broadPeak  -f 0.5 -wa > ${PeakDir}/Input_consensus.narrowPeak
bedtools intersect -a ${PeakDir}/${dpf3_H3K27me3_1}_peaks.broadPeak -b ${PeakDir}/${dpf3_H3K27me3_2}_peaks.broadPeak ${PeakDir}/${dpf3_H3K27me3_3}_broadPeak.bed -f 0.5 -wa > ${PeakDir}/dpf3_H3K27me3_Consensus.broadPeak
bedtools intersect -a ${PeakDir}/${dpf6_H3K27me3_1}_peaks.broadPeak -b ${PeakDir}/${dpf6_H3K27me3_2}_peaks.broadPeak ${PeakDir}/${dpf6_H3K27me3_3}_broadPeak.bed -f 0.5 -wa > ${PeakDir}/dpf6_H3K27me3_Consensus.broadPeak
bedtools intersect -a ${PeakDir}/${Mycelia_H3K27me3_1}_peaks.broadPeak -b ${PeakDir}/${Mycelia_H3K27me3_2}_peaks.broadPeak -f 0.5 -wa > ${PeakDir}/Mycelia_H3K27me3_Consensus.broadPeak

## Intersect consensus peakset (made in R) with N. crassa promoters
# First, make a bed file adding 500bp to the start site of each gene
bedtools slop -s  -i ~/NcGenome/NcGenes.bed -g ~/NcGenome/chrom_sizes.txt -l 500 > ~/NcGenome/NcGenes_StartMinus500.bed
bedtools intersect -a ~/NcGenome/NcGenes_StartMinus500.bed -b all_H3K27me3_Peaks.bed -f 0.25 -wa > H3K27me3_Gene_Promoters.bed





## report no. peaks in each peakfile
wc -l * > peak_counts.txt

## for each peakfile, calculate peak coverage
GENOME_SIZE=38639769

for f in *; do
    awk -v file="$f" -v gsize="$GENOME_SIZE" '
        $3 > $2 { sum += ($3 - $2) }
        END {
            pct = (sum / gsize) * 100
            printf "%s\t%d\t%.4f%%\n", file, sum, pct
        }
    ' "$f"
done


## remove background peaks
bedtools intersect -a  ${PeakDir}/fsd1_Consensus_Peaks_raw.narrowPeak -b ${PeakDir}/GFPtrap_Consensus_Peaks.narrowPeak -f 0.8 -v >  ${PeakDir}/fsd1_Consensus_Peaks.bed
bedtools intersect -a  ${PeakDir}/h2aZ_M_Consensus_Peaks_raw.broadPeak -b ${PeakDir}/GFPtrap_Consensus_Peaks.narrowPeak -f 0.8 -v >  ${PeakDir}/h2aZ_M_Consensus_Peaks.bed
bedtools intersect -a  ${PeakDir}/h2aZ_dpf3_Consensus_Peaks_raw.broadPeak -b ${PeakDir}/GFPtrap_Consensus_Peaks.narrowPeak -f 0.8 -v >  ${PeakDir}/h2aZ_dpf3_Consensus_Peaks.bed
bedtools intersect -a  ${PeakDir}/h2aZ_dpf6_Consensus_Peaks_raw.broadPeak -b ${PeakDir}/GFPtrap_Consensus_Peaks.narrowPeak -f 0.8 -v >  ${PeakDir}/h2aZ_dpf6_Consensus_Peaks.bed

## Use homer to make a consensus peakset
 ml Homer/5.1-foss-2023a-R-4.3.2

mergePeaks -d 100 ${PeakDir}/${Control1}_peaks.narrowPeak ${PeakDir}/${Control2}_peaks.narrowPeak ${PeakDir}/${Control3}_peaks.narrowPeak -venn ${PeakDir}/control_peaks.txt  > ${PeakDir}/WT_GFPtrap_Peaks.bed
mergePeaks -d 100 ${PeakDir}/${Chip1}_peaks.narrowPeak ${PeakDir}/${Chip2}_peaks.narrowPeak -venn ${PeakDir}/fsd1_peaks.txt -prefix fsd1 > ${PeakDir}/fsd1_GFPtrap_Peaks.bed

findMotifsGenome.pl ${PeakDir}/fsd1_GFPtrap_Peaks.bed /home/ad45368/NcGenome/GCA_000182925.2_NC12_genomic.fna ${Motifs}/ -size given -bg ${PeakDir}/WT_GFPtrap_Peaks.bed

ml ucsc/443

bigWigMerge  ${bwdir}/151-130_ChIP_WT_M_H3K27me3_Rep2_S153_L002.bin_25.smooth_75Bulk.bw ${bwdir}/150-31_ChIP_WT_M_H3K27me3_Rep1_S31_L002.bin_25.smooth_75Bulk.bw  ${bwdir}/MyceliaH3K27me3_merge.bedGraph
bedGraphToBigWig ${bwdir}/MyceliaH3K27me3_merge.bedGraph /home/ad45368/chrom_sizes.txt  ${bwdir}/MyceliaH3K27me3_merge.bw

bigWigMerge  ${bwdir}/150-28_ChIP_WT_P_6pf_H3K27me3_Rep2_S28_L002.bin_25.smooth_75Bulk.bw ${bwdir}/152-38_ChIP_WT_6dpf_P_H3K27me3_Rep5_S38_L006.bin_25.smooth_75Bulk.bw ${bwdir}/dpf6H3K27me3_merge.bedGraph
bedGraphToBigWig ${bwdir}/dpf6H3K27me3_merge.bedGraph /home/ad45368/chrom_sizes.txt  ${bwdir}/dpf6H3K27me3_merge.bw

ml ucsc/443

bigWigMerge ${bwdir}/150-26_ChIP_WT_P_3dpf_H3K36me3_Rep1_S26_L002.bin_25.smooth_75Bulk.bw ${bwdir}/151-127_ChIP_WT_P_3dpf_H3K36me3_Rep2_.bin_25.smooth_75Bulk.bw ${bwdir}/dpf3H3K36me3_merge.bedGraph
bigWigMerge ${bwdir}/150-32_ChIP_WT_M_H3K36me3_Rep1_S32_L002.bin_25.smooth_75Bulk.bw ${bwdir}/152-30_ChIP_WT_H3K36me3_Rep4_S30_L006.bin_25.smooth_75Bulk.bw ${bwdir}/MyceliaH3K36me3_merge.bedGraph
bigWigMerge ${bwdir}/152-39_ChIP_WT_6dpf_P_H3K36me3_Rep5_S39_L006.bin_25.smooth_75Bulk.bw ${bwdir}/151-129_ChIP_WT_P_6pf_H3K36me3_Rep3_.bin_25.smooth_75Bulk.bw ${bwdir}/dpf6H3K36me3_merge.bedGraph

bedGraphToBigWig ${bwdir}/dpf3H3K36me3_merge.bedGraph /home/ad45368/chrom_sizes.txt  ${bwdir}/dpf3H3K36me3_merge.bw
bedGraphToBigWig ${bwdir}/dpf6H3K36me3_merge.bedGraph /home/ad45368/chrom_sizes.txt  ${bwdir}/dpf6H3K36me3_merge.bw
bedGraphToBigWig ${bwdir}/MyceliaH3K36me3_merge.bedGraph /home/ad45368/chrom_sizes.txt  ${bwdir}/MyceliaH3K36me3_merge.bw

ml deepTools/3.5.5-gfbf-2023a
base="/lustre2/scratch/ad45368/ChIPseq/FinalH3k27me3_ChIP"
gfp_base="/lustre2/scratch/ad45368/ChIPseq/GFP_ChIP/MappingOutput/bigWig"
meta=${base}/Meta
ml deepTools/3.5.5-gfbf-2023a

computeMatrix reference-point --referencePoint TSS  -R ${base}/SharedK27Genes.bed ${base}/PeritheciaOnlyK27Genes.bed ${base}/MyceliaOnlyK27Genes.bed -b 500 -a 300   \
    -S ${base}/bigWig/MyceliaH3K27me3_merge.bw ${base}/bigWig_H3K36/MyceliaH3K36me3_merge.bw  ${base}/bigWig/SRR5177530.bw ${base}/bigWig/dpf6H3K27me3_merge.bw ${base}/bigWig_H3K36/dpf3H3K36me3_merge.bw ${base}/bigWig_H3K36/dpf6H3K36me3_merge.bw ${base}/bigWig/SRR5177522.bw \
    --samplesLabel "Mycelia H3K27me3" "Mycelia H3K36me3"  "Mycelia RNA" "dpf6 H3K27me3" "dpf3 H3K36me3" "dpf6 H3K36me3" "dpf6 RNA" \
    -o ${meta}/h3k27genes.gz \
    --outFileSortedRegions ${meta}/h3k27me3_genes_mat.bed \
    -p 6

computeMatrix reference-point --referencePoint TSS  -R ${base}/SharedK27Genes.bed ${base}/PeritheciaOnlyK27Genes.bed ${base}/MyceliaOnlyK27Genes.bed -b 500 -a 300   \
    -S ${base}/bigWig/MyceliaH3K27me3_merge.bw ${base}/bigWig_H3K36/MyceliaH3K36me3_merge.bw /lustre2/scratch/ad45368/ChIPseq/GFP_ChIP/MappingOutput/bigWig/153_51_ChIP_h2aZ_M_GFP_Rep1_.bin_25.smooth_75Bulk.bw  ${base}/bigWig/WT_mycelia_log.bw ${base}/bigWig/dpf6H3K27me3_merge.bw /lustre2/scratch/ad45368/ChIPseq/GFP_ChIP/MappingOutput/bigWig/153_57_ChIP_h2aZ_dpf6_GFP_Rep1_.bin_25.smooth_75Bulk.bw ${base}/bigWig_H3K36/dpf6H3K36me3_merge.bw ${base}/bigWig/SRR5177521_log.bw \
    --samplesLabel "Mycelia H3K27me3" "Mycelia H3K36me3" "Mycelia H2A.Z" "Mycelia RNA" "dpf6 H3K27me3" "dpf6 H3K36me3" "dpf6 H2A.Z" "dpf6 RNA" \
    -o ${meta}/h3k27genes_H2AZ.gz \
    --missingDataAsZero \
    --outFileSortedRegions ${meta}/h3k27me3_genes_H2AZ_mat.bed \
    -p 6

   plotHeatmap --matrixFile ${meta}/h3k27genes_H2AZ.bed -o ${meta}/h3k27me3_genes_H2AZ_mat.png --outFileNameMatrix ${meta}/h3k27me3_genes_H2AZ_matout.gz       --colorMap Greens Oranges Blues Reds Greens Oranges Blues Reds       --zMax 70 30 10 60 70 30 10   --missingDataColor white --heatmapHeight 20

computeMatrix reference-point --referencePoint TSS  -R ${base}/SharedK27Genes.bed ${base}/PeritheciaOnlyK27Genes.bed ${base}/MyceliaOnlyK27Genes.bed -b 500 -a 300   \
    -S ${base}/bigWig/MyceliaH3K27me3_merge.bw ${base}/bigWig_H3K36/MyceliaH3K36me3_merge.bw  ${base}/bigWig/WT_mycelia_log.bw ${base}/bigWig/dpf6H3K27me3_merge.bw ${base}/bigWig_H3K36/dpf3H3K36me3_merge.bw ${base}/bigWig_H3K36/dpf6H3K36me3_merge.bw ${base}/bigWig/SRR5177521_log.bw \
    --samplesLabel "Mycelia H3K27me3" "Mycelia H3K36me3"  "Mycelia RNA" "dpf6 H3K27me3" "dpf3 H3K36me3" "dpf6 H3K36me3" "dpf6 RNA" \
    -o ${meta}/h3k27genes.gz \
    --outFileSortedRegions ${meta}/h3k27me3_genes_mat.bed \
    -p 6

   plotHeatmap --matrixFile ${meta}/h3k27genes.gz -o ${meta}/h3k27me3_gene.png --outFileNameMatrix ${meta}/h3k27me3_genes_matout.gz       --colorMap Greens Oranges Reds Greens Oranges Oranges Reds       --zMax 70 30 10 60 70 30 10   --missingDataColor white --heatmapHeight 20


