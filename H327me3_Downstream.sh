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

bigWigMerge  ${bwdir}/${Control1}.bin_25.smooth_75Bulk.bw ${bwdir}/${Control2}.bin_25.smooth_75Bulk.bw ${bwdir}/${Control3}.bin_25.smooth_75Bulk.bw ${bwdir}/GFPtrap_Control_merge.bedGraph
bedGraphToBigWig ${bwdir}/GFPtrap_Control_merge.bedGraph /home/ad45368/chrom_sizes.txt  ${bwdir}/GFPtrap_Control_merge.bw

bigWigMerge  ${bwdir}/${Chip1}.bin_25.smooth_75Bulk.bw ${bwdir}/${Chip2}.bin_25.smooth_75Bulk.bw ${bwdir}/fsd1_dpf6_merged.bedGraph
bedGraphToBigWig ${bwdir}/fsd1_dpf6_merged.bedGraph /home/ad45368/chrom_sizes.txt  ${bwdir}/fsd1_dpf6_merged.bw
