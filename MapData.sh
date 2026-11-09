#!/bin/bash
# MapData.sh  –  CUT&RUN mapping pipeline
# Adapted for UW-Madison CHTC / HTCondor + Apptainer containers
#
# Variables injected by HTCondor via the submit file:
#   accession   – sample name (e.g. SRR12345678)
#   fastqPath   – path to FastqFiles directory (accessible on execute node)
#   outdir      – path to MappingOutput directory
#   GENOME      – path to BWA-indexed genome (must be staged or on OSDF)
#   THREADS     – number of CPU threads  (matches request_cpus in submit file)
#   BIN         – deepTools bin size
#   SMOOTH      – deepTools smooth length
set -euo pipefail
 
echo "===== Starting CUT&RUN mapping for: ${accession} ====="
echo "Host: $(hostname)"
echo "Date: $(date)"

###################
#start
####################################################

#input file variables
  read1=(${fastqPath}/${accession}*_R1_001.fastq.gz)
  read2=(${fastqPath}/${accession}*_R2_001.fastq.gz)

   
if [[ ! -f "${read1[0]}" ]]; then
    echo "ERROR: R1 fastq not found for ${accession} in ${fastqPath}" >&2
    exit 1
fi
if [[ ! -f "${read2[0]}" ]]; then
    echo "ERROR: R2 fastq not found for ${accession} in ${fastqPath}" >&2
    exit 1
fi
###################################


#make output file folders
trimmed="${outdir}/TrimmedFastQs/${accession}"
mkdir $trimmed

tmp="${outdir}/tempFile"
mkdir -p "${tmp}"

bamdir="${outdir}/bamFiles"
mkdir -p "${bamdir}"

bwDir="${outdir}/bigWig"
mkdir "${bwDir}"

name=$(echo "$accession" | sed -E 's/_S[0-9]{1,3}_L[0-9]{3}//')

PeakDir="${outdir}/Peaks"
mkdir ${PeakDir}

#make variables for output file names
bam="${bamdir}/${name}.bam"
bigwig="${bwDir}/${name}"
peak="$PeakDir/${name}"


############# Read Trimming ##############
#remove adaptors, trim low quality reads (default = phred 20), length > 25

##fastq files from the ebi link are in folders that either have one file with a SRR##.fastq.gz or a SRR##_1.fastq.gz ending, or have two files with a SRR##_1.fastq.gz ending or a SRR##_2.fastq.gz ending
#or have three files with a SRR##_1.fastq.gz, SRR##_2.fastq.gz and SRR##.fastq.gz ending. In this case, the third file corresponds to unpaired reads that the depositers mapped.

#JGI files should have a SRR##_1.fastq.gz ending and a SRR##_2.fastq.gz file

##################
#Trimming
#################
	  module load Trim_Galore/0.6.10-GCCcore-12.3.0

	  trim_galore --illumina --fastqc --paired --length 25 --basename ${accession} --gzip -o $trimmed $read1 $read2
	  wait


module load  SAMtools/1.21-GCC-13.3.0
module load  BWA/0.7.18-GCCcore-13.3.0

#

#make directory to store temporary files written by samtools sort
mkdir -p ${tmp}/${accession}

bwa mem -M -v 3 -t $THREADS $GENOME ${trimmed}/*val_1.fq.gz ${trimmed}/*val_2.fq.gz | samtools view -bhSu - | samtools sort -@ $THREADS -T ${tmp}/${accession} -o "$bam" -
samtools index "$bam"

#delete directory written by samtools sort
rm -r ${tmp}/${accession}

module load deepTools/3.5.5-gfbf-2023a
#Plot all reads
bamCoverage -p $THREADS -bs $BIN --normalizeUsing BPM --minMappingQuality 20 --smoothLength $SMOOTH -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}Bulk.bw"

#plot mononucleosomes
bamCoverage -p $THREADS --MNase -bs 1 --normalizeUsing BPM --minMappingQuality 20 --smoothLength 25 -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}_MNase.bw"

#call Peaks
module load MACS3/3.0.1-gfbf-2023a

#using --nolambda paramenter to call peaks without control
macs3 callpeak -t "${bam}" -f BAMPE -n "${accession}" --broad -g 41037538 --broad-cutoff 0.1 --outdir "${PeakDir}" --min-length 800 --max-gap 500 --nolambda

echo "===== Finished: ${accession} ====="



