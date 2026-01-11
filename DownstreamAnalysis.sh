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
outdir="../MappingOutput"
bamdir="${outdir}/bamFiles"
bwdir="${outdir}/bigWig"
PeakDir="${outdir}/Peaks"
Motfis="${outdir}/Motifs"
## Input control and rep information, using name:

Control1="143-144_ChIP_WT_gfp-trap_Rep1"
Control2="153-44_ChIP_WT_dpf3_gfp-trap_Rep1"
Control3="153-47_ChIP_WT_dpf6_gfp-trap_Rep1"
Chip1="143-3_ChIP_fsd1-gfpXsad1_GFPtrap_Rep1"
Chip2="143-5_ChIP_fsd1-gfpXsad1_GFPtrap_Rep1"

KC683708.1
## Use homer to make a consensus peakset
 ml Homer/5.1-foss-2023a-R-4.3.2

mergePeaks -d 100 ${PeakDir}/${Control1}_peaks.narrowPeak ${PeakDir}/${Control2}_peaks.narrowPeak ${PeakDir}/${Control3}_peaks.narrowPeak -venn ${PeakDir}/control_peaks.txt  > ${PeakDir}/WT_GFPtrap_Peaks.bed
mergePeaks -d 100 ${PeakDir}/${Chip1}_peaks.narrowPeak ${PeakDir}/${Chip2}_peaks.narrowPeak -venn ${PeakDir}/fsd1_peaks.txt -prefix fsd1 > ${PeakDir}/fsd1_GFPtrap_Peaks.bed

findMotifsGenome.pl ${PeakDir}/fsd1_GFPtrap_Peaks.bed /home/ad45368/NcGenome/GCA_000182925.2_NC12_genomic.fna ${Motifs}/ -size given -bg ${PeakDir}/WT_GFPtrap_Peaks.bed

ml ucsc/443

bigWigMerge  ${bwdir}/${Control1}.bin_25.smooth_75Bulk.bw ${bwdir}/${Control2}.bin_25.smooth_75Bulk.bw ${bwdir}/${Control3}.bin_25.smooth_75Bulk.bw ${bwdir}/GFPtrap_Control_merge.bedGraph
bedGraphToBigWig ${bwdir}/GFPtrap_Control_merge.bedGraph /home/ad45368/chrom_sizes.txt  ${bwdir}/GFPtrap_Control_merge.bw

bigWigMerge  ${bwdir}/${Chip1}.bin_25.smooth_75Bulk.bw ${bwdir}/${Chip2}.bin_25.smooth_75Bulk.bw ${bwdir}/fsd1_dpf6_merged.bedGraph
bedGraphToBigWig ${bwdir}/fsd1_dpf6_merged.bedGraph /home/ad45368/chrom_sizes.txt  ${bwdir}/fsd1_dpf6_merged.bw
