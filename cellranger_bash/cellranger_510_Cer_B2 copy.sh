#!/bin/bash
#SBATCH --job-name=cellranger_510_Cer_B2
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --ntasks=8
#SBATCH --output=cellranger_510_Cer_B2_%j.out
#SBATCH --error=cellranger_510_Cer_B2_%j.err

module purge
module load bluebear
module load bear-apps/2023a
module load CellRanger/9.0.1
set -e

cellranger count --id=510_Cer_B2 \
  --transcriptome=/rds/projects/c/colemdjl-exphaem/MSc_DL_Murali/rawdata/refdata-gex-GRCh38-2024-A \
  --fastqs=/rds/projects/c/colemdjl-exphaem/MSc_DL_Murali/rawdata/510_Cer/B2 \
  --sample=510_Cer_B2_RNA \
  --localcores=8 \
  --localmem=32 \
  --create-bam false
