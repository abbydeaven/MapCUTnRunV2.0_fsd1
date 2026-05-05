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
outdir="/scratch/ad45368/ChIPseq/GFP_ChIP/MappingOutput"
bamdir="${outdir}/bamFiles"
bwdir="${outdir}/bigWig"
PeakDir="${outdir}/Peaks"
Motifs="${outdir}/Motifs"
meta="${outdir}/Metaplots"
## Input control and rep information, using name:

#Control1="143_144_ChIP_WT_gfp_trap_Rep1"
Control2="153_44_ChIP_WT_dpf3_gfp_trap_Rep1"
Control3="153_47_ChIP_WT_dpf6_gfp_trap_Rep1"
fsd1_1="143_3_ChIP_fsd1_gfpXsad1_GFPtrap_Rep1"
fsd1_2="143_5_ChIP_fsd1_gfpXsad1_GFPtrap_Rep1"

h2aZ_gt_M1="153_50_ChIP_h2aZ_M_gfp_trap_Rep1"
h2aZ_gp_M1="153_51_ChIP_h2aZ_M_GFP_Rep1"

h2aZ_gp_d3="153_54_ChIP_h2aZ_dpf3_GFP_Rep1"
h2aZ_gt_d3="153_53_ChIP_h2aZ_dpf3_gfp_trap_Rep1"

h2aZ_gt_d6="153_56_ChIP_h2aZ_dpf6_gfp_trap_Rep1"
h2aZ_gp_d6="153_57_ChIP_h2aZ_dpf6_GFP_Rep1"

for f in *; do
      name=$(echo "$f" | sed -E 's/\-/\_/')
      mv $f $name
      done

## Use bedtools to make a consensus peakset of peaks with 80% overlap in all samples
ml BEDTools/2.31.1-GCC-13.3.0

bedtools intersect -a $PeakDir/${Control2}.broadPeak -b ${PeakDir}/${Control3}.broadPeak  -f 0.5 -wa > ${PeakDir}/GFPtrap_Consensus_Peaks.narrowPeak
bedtools intersect -a ${PeakDir}/${fsd1_1}.narrowPeak -b ${PeakDir}/${fsd1_1}.narrowPeak -f 0.5 -wa > ${PeakDir}/fsd1_Consensus_Peaks.bed
bedtools intersect -a ${PeakDir}/${h2aZ_gt_M1}.broadPeak -b ${PeakDir}/${h2aZ_gp_M1}.broadPeak -f 0.5 -wa > ${PeakDir}/h2aZ_M_Consensus_Peaks_raw.broadPeak
bedtools intersect -a ${PeakDir}/${h2aZ_gp_d3}.broadPeak -b ${PeakDir}/${h2aZ_gt_d3}.broadPeak -f 0.5 -wa > ${PeakDir}/h2aZ_dpf3_Consensus_Peaks_raw.broadPeak
bedtools intersect -a ${PeakDir}/${h2aZ_gt_d6}.broadPeak -b ${PeakDir}/${h2aZ_gp_d6}.broadPeak -f 0.5 -wa > ${PeakDir}/h2aZ_dpf6_Consensus_Peaks_raw.broadPeak

## report no. peaks in each peakfile
wc -l * > peak_counts.txt


## for each peakfile, calculate peak coverage
GENOME_SIZE=41037538 

for f in *; do
    awk -v file="$f" -v gsize="$GENOME_SIZE" '
        $3 > $2 { sum += ($3 - $2) }
        END {
            pct = (sum / gsize) * 100
            printf "%s\t%d\t%.4f%%\n", file, sum, pct
        }
    ' "$f"
done


## remove background peaks. might not have to do this.
#bedtools intersect -a  ${PeakDir}/fsd1_Consensus_Peaks_raw.narrowPeak -b ${PeakDir}/GFPtrap_Consensus_Peaks.narrowPeak -f 0.9 -v >  ${PeakDir}/fsd1_Consensus_Peaks.bed
bedtools intersect -a  ${PeakDir}/h2aZ_M_Consensus_Peaks_raw.broadPeak -b ${PeakDir}/GFPtrap_Consensus_Peaks.narrowPeak -f 0.9 -v >  ${PeakDir}/h2aZ_M_Consensus_Peaks.bed
bedtools intersect -a  ${PeakDir}/h2aZ_dpf3_Consensus_Peaks_raw.broadPeak -b ${PeakDir}/GFPtrap_Consensus_Peaks.narrowPeak -f 0.8 -v >  ${PeakDir}/h2aZ_dpf3_Consensus_Peaks.bed
bedtools intersect -a  ${PeakDir}/h2aZ_dpf6_Consensus_Peaks_raw.broadPeak -b ${PeakDir}/GFPtrap_Consensus_Peaks.narrowPeak -f 0.8 -v >  ${PeakDir}/h2aZ_dpf6_Consensus_Peaks.bed

## Use homer to identify motifs (findMotifsGenome). then, run homer to identify OTHER TF motifs in fsd-1 target promoters.
 ml Homer/5.1-foss-2023a-R-4.3.2

#mergePeaks -d 100 ${PeakDir}/${Control1}_peaks.narrowPeak ${PeakDir}/${Control2}_peaks.narrowPeak ${PeakDir}/${Control3}_peaks.narrowPeak -venn ${PeakDir}/control_peaks.txt  > ${PeakDir}/WT_GFPtrap_Peaks.bed
#mergePeaks -d 100 ${PeakDir}/${Chip1}_peaks.narrowPeak ${PeakDir}/${Chip2}_peaks.narrowPeak -venn ${PeakDir}/fsd1_peaks.txt -prefix fsd1 > ${PeakDir}/fsd1_GFPtrap_Peaks.bed

## IMPORTANT!! downgrade to an earlier version of perl for compatibility.
ml Perl/5.34.0-GCCcore-11.2.0
findMotifsGenome.pl ${PeakDir}/fsd1_peaks_devup.bed Ncrassa ${Motifs}/Development_Up -size 50 -bg ${PeakDir}/GFPtrap_Consensus_Peaks.narrowPeak
findMotifsGenome.pl ${PeakDir}/fsd1_peaks_devdown.bed Ncrassa ${Motifs}/Development_Down -size 50 -bg ${PeakDir}/GFPtrap_Consensus_Peaks.narrowPeak
findMotifsGenome.pl ${PeakDir}/fsd1_peaks_nondev.bed Ncrassa ${Motifs}/NonDevelopment -size 50 -bg ${PeakDir}/GFPtrap_Consensus_Peaks.narrowPeak

# annotate peaks + identify where motifs are. motif1 = highest prevalence (~10%), most similar to ndt80
/scratch/ad45368/ChIPseq/GFP_ChIP/MappingOutput/Peaks/fsd1_Consensus_Peaks.bed Ncrassa /scratch/ad45368/ChIPseq/GFP_ChIP/MappingOutput/Motifs/ -size 200 -bg /scratch/ad45368/ChIPseq/GFP_ChIP/MappingOutput/Peaks/GFPtrap_Consensus_Peaks.narrowPeak
annotatePeaks.pl ${PeakDir}/fsd1_Consensus_Peaks.bed Ncrassa -gtf ~/NcGenome/fungiDB_GFFtoGTF_conversion.gtf > ${Motifs}/fsd1_annotated_peaks.txt

ml ucsc/443

bigWigMerge  ${bwdir}/153_44_ChIP_WT_dpf3_gfp_trap_Rep1_S44_L002.bin_25.smooth_75Bulk.bw ${bwdir}/${Control3}.bin_25.smooth_75Bulk.bw ${bwdir}/GFPtrap_Control_merge.bedGraph
bedGraphToBigWig ${bwdir}/GFPtrap_Control_merge.bedGraph /home/ad45368/chrom_sizes.txt  ${bwdir}/GFPtrap_Control_merge.bw

bigWigMerge  ${bwdir}/${fsd1_1}.bin_25.smooth_75Bulk.bw ${bwdir}/${fsd1_2}.bin_25.smooth_75Bulk.bw ${bwdir}/fsd1_dpf6_merged.bedGraph
bedGraphToBigWig ${bwdir}/fsd1_dpf6_merged.bedGraph /home/ad45368/chrom_sizes.txt  ${bwdir}/fsd1_dpf6_merged.bw


## make fsd1 sashimi plots
. env/ggsashimi/bin/activate
 ~/ggsashimi/ggsashimi.py -b fsd1_dev.tsv -c CM002239.1:5296170-5303095 -o fsd1_dev_2 -S plus -M 5 -g ~/NcGenome/fungiDB_GFFtoGTF_conversion.gtf -O 3 -C 3 -L 1 --fix-y-scale --ann-height 1 --base-size 18 -F png
[W::hts_idx_load3] The index file is older than the data file: ../MappingOutput/Temp/SRR5177527_tmp3.bam.bai

