#!/bin/bash
#SBATCH --partition=batch
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ad45368@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100gb
#SBATCH --time=8:00:00
#SBATCH --output=../ChIPOutput/logs/%x.out
#SBATCH --error=../ChIPOutput/logs/%x.err


cd $SLURM_SUBMIT_DIR
##Directory information+variables:
outdir="../ChIPOutput"
bam="${outdir}/bamFiles"
bw="${outdir}/bigWig"
PeakDir="${outdir}/Peaks"
Motifs="${outdir}/Motifs"
meta="${outdir}/Metaplots"
## Input control and rep information, using name:

#Control1="143_144_ChIP_WT_gfp_trap_Rep1"
Control1="143-144_ChIP_WT_gfp-trap_Rep1"
Control2="153-47_ChIP_WT_dpf6_gfp-trap_Rep1"
fsd1_1="143_3_ChIP_fsd1_gfpXsad1_GFPtrap_Rep1"
fsd1_2="143_5_ChIP_fsd1_gfpXsad1_GFPtrap_Rep1"

# Replace '-' in filenames with "_" to avoid issues with downstream analysis
for f in *; do
      name=$(echo "$f" | sed -E 's/\-/\_/')
      mv $f $name
      done


## Use bedtools to make a consensus peakset of peaks with 80% overlap in all samples
ml BEDTools/2.31.1-GCC-13.3.0

bedtools intersect -a $PeakDir/${Control1}.narrowPeak -b ${PeakDir}/${Control2}.narrowPeak  -f 0.5 -wa > ${PeakDir}/GFPtrap_Consensus_Peaks.bed
bedtools intersect -a ${PeakDir}/${fsd1_1}.narrowPeak -b ${PeakDir}/${fsd1_1}.narrowPeak -f 0.5 -wa > ${PeakDir}/fsd1_Consensus_Peaks.bed

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

## merge bigWig files 
ml ucsc/443

bigWigMerge  ${bwdir}/${Control1}.bin_25.smooth_75Bulk.bw ${bwdir}/${Control2}.bin_25.smooth_75Bulk.bw ${bwdir}/GFPtrap_Control_merge.bedGraph
bedGraphToBigWig ${bwdir}/GFPtrap_Control_merge.bedGraph /home/ad45368/chrom_sizes.txt ${bwdir}/GFPtrap_Control_merge.bw

bigWigMerge  ${bwdir}/${fsd1_1}.bin_25.smooth_75Bulk.bw ${bwdir}/${fsd1_2}.bin_25.smooth_75Bulk.bw ${bwdir}/fsd1_dpf6_merged.bedGraph
bedGraphToBigWig ${bwdir}/fsd1_dpf6_merged.bedGraph /home/ad45368/chrom_sizes.txt  ${bwdir}/fsd1_dpf6_merged.bw


## Identify motifs using MEME suite
ml MEME/5.5.7-gompi-2023b

    # First convert peak files to FASTA using bed2fasta
    bed2fasta ${PeakDir}/GFPtrap_Consensus_Peaks.bed /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic.fna -o ${Motifs}/GFPtrap_background_peaksequences.fasta
    bed2fasta ${PeakDir}/fsd1_Consensus_Peaks.bed /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic.fna -o ${Motifs}/fsd1_peaksequences.fasta

    # Identify motifs
    meme-chip ${Motifs}/fsd1_peaksequences.fasta \
    -neg ${Motifs}/GFPtrap_background_peaksequences.fasta \
    -oc ${Motifs}/fsd1_motif_analysis \
    -meme-mod zoops \
    -meme-nmotifs 5
