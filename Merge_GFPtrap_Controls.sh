#!/bin/bash
#SBATCH --job-name=merge_chip_control
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32gb
#SBATCH --time=2:00:00
#SBATCH --output=../ChIPOutput/logs/merge_control.%j.out
#SBATCH --error=../ChIPOutput/logs/merge_control.%j.err


cd "${SLURM_SUBMIT_DIR}"

## Make variable names and directories
outdir="../ChIPOutput"
control_name="GFPtrap_control"

bamdir="${outdir}/bamFiles"

module load SAMtools/1.21-GCC-13.3.0

samtools merge -@ 4 -o $bamdir/GFPtrap_control.bam "${bamdir}/143-144_ChIP_WT_gfp-trap_Rep1.bam" "${bamdir}/153-47_ChIP_WT_dpf6_gfp-trap_Rep1.bam"
samtools index -@ 4 $bamdir/GFPtrap_control.bam
