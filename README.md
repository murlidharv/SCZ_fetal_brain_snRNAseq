# SCZ_fetal_brain_snRNAseq

Reanalysis of human fetal brain single-nucleus RNA sequencing (snRNA-seq) data examining region- and cell-type-specific expression of convergent schizophrenia (SCZ) risk genes during mid-gestational neurodevelopment.

## Study Overview

This project investigates the developmental expression patterns of four convergent schizophrenia-associated genes:

- SP4
- GRIN2A
- STAG1
- KLC1

across multiple fetal brain regions, with particular focus on radial glial (RG/oRG) and neuronal populations.

## Dataset

Publicly available human fetal brain snRNA-seq data:

- Cameron et al., 2023
- European Genome-Phenome Archive (EGA)
- Accession: EGAD00001009303

Brain regions analysed:

- Frontal Cortex (FC)
- Ganglionic Eminence (GE)
- Hippocampus (Hipp)
- Thalamus (Thal)
- Cerebellum (Cer)

## Main Analysis Pipeline

1. Quality control and preprocessing
2. Seurat-based integration
3. Clustering and cell-type annotation
4. oRG reannotation using canonical markers
5. Differential expression analysis
6. Final figure generation and visualisation

## Repository Structure

- `scripts/` — R analysis scripts
- `cellranger_bash/` — representative SLURM Cell Ranger scripts
- `fastq_metadata/` — FASTQ summary files and metadata
- `marker_gene_lists/` — top cluster marker gene tables
- `docs/` — workflow figures and supporting documentation

## Software

- R
- Seurat v5
- Cell Ranger 9.0.1
- scDblFinder

## Differential Expression

Differential expression between oRG and neuronal populations was performed using Seurat `FindMarkers` with Wilcoxon rank-sum testing and FDR correction.

## Notes

This repository contains selected scripts and supporting documentation only.

Raw FASTQ files, Seurat objects, controlled-access EGA data, and credential files are not included in this public repository.

## Author

Venkiteswaran Muralidhar  
University of Birmingham  
MSc Bioinformatics
