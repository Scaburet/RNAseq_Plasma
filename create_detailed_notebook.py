import nbformat
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
import json

def create_genome_indexing_cells():
    """Create cells for genome downloading and indexing section"""
    cells = []

    # Add documentation
    doc_md = """## Reference Genome Download and Indexing

### Important Notes
- The mouse reference genome (GRCm39) will be downloaded from Ensembl
- We'll use STAR for indexing the genome
- This step requires significant computational resources
- The index will be used for all subsequent mapping steps

### Key Parameters for STAR Indexing
- `--runMode genomeGenerate`: Tells STAR to create an index
- `--genomeDir`: Directory where the index will be stored
- `--genomeFastaFiles`: Path to the reference genome FASTA
- `--sjdbGTFfile`: Path to the GTF annotation file
- `--runThreadN`: Number of threads (we'll use 10)

### Documentation Links
- [Ensembl Mouse Genome](https://www.ensembl.org/Mus_musculus/Info/Index)
- [STAR Manual - Genome Generation](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf#page=7)
"""
    cells.append(new_markdown_cell(doc_md))

    # Add genome download code
    download_code = """# Download and prepare reference genome
mkdir -p reference
cd reference

# Download mouse reference genome from Ensembl
wget https://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz

# Create STAR index
mkdir -p star_index
STAR --runMode genomeGenerate \\
     --genomeDir star_index \\
     --genomeFastaFiles Mus_musculus.GRCm39.dna.primary_assembly.fa \\
     --sjdbGTFfile /srv/data/meg-m2-rnaseq/Genomes/Mmu/GRCm39/extracted/mouse\ genome_annotation-M35.gtf \\
     --runThreadN 10"""
    cells.append(new_code_cell(download_code))
    return cells

def create_quality_check_cells():
    """Create cells for quality check section"""
    cells = []

    # Add documentation
    doc_md = """## Quality Control with FastQC

### Overview
FastQC performs quality control checks on raw sequence data. We'll analyze:
- Basic Statistics
- Per base sequence quality
- Per sequence quality scores
- Sequence duplication levels
- Overrepresented sequences

### Important Parameters
- `-o`: Output directory
- `-t`: Number of threads
- `--noextract`: Don't extract zip files

### Documentation
- [FastQC Documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Example Reports](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html)
"""
    cells.append(new_markdown_cell(doc_md))

    # Add FastQC code
    fastqc_code = """# Run FastQC on first two samples
DATA_DIR="/srv/data/meg-m2-rnaseq/Data/fastq/raw"
fastqc -o Results/fastqc -t 10 \\
    ${DATA_DIR}/sample1_R1.fastq.gz \\
    ${DATA_DIR}/sample1_R2.fastq.gz \\
    ${DATA_DIR}/sample2_R1.fastq.gz \\
    ${DATA_DIR}/sample2_R2.fastq.gz

# Command for all samples (shown but not executed)
: '
fastqc -o Results/fastqc -t 10 ${DATA_DIR}/*.fastq.gz
'"""
    cells.append(new_code_cell(fastqc_code))
    return cells

def create_preprocessing_cells():
    """Create cells for preprocessing section"""
    cells = []

    # Add documentation
    doc_md = """## Read Preprocessing with fastp

### Overview
fastp performs quality control and preprocessing:
- Adapter trimming
- Quality filtering
- Base correction
- Quality reporting

### Key Parameters
- `--in1/--in2`: Input paired-end files
- `--out1/--out2`: Output paired-end files
- `--html`: HTML report output
- `--json`: JSON report output
- `--thread`: Number of threads

### Documentation
- [fastp GitHub](https://github.com/OpenGene/fastp)
- [Parameters Guide](https://github.com/OpenGene/fastp#all-options)
"""
    cells.append(new_markdown_cell(doc_md))

    # Add fastp code
    fastp_code = """# Process first two samples with fastp
DATA_DIR="/srv/data/meg-m2-rnaseq/Data/fastq/raw"
for sample in sample1 sample2; do
    fastp --in1 ${DATA_DIR}/${sample}_R1.fastq.gz \\
          --in2 ${DATA_DIR}/${sample}_R2.fastq.gz \\
          --out1 Results/fastp/${sample}_R1_cleaned.fastq.gz \\
          --out2 Results/fastp/${sample}_R2_cleaned.fastq.gz \\
          --html Results/fastp/${sample}_report.html \\
          --json Results/fastp/${sample}_report.json \\
          --thread 10
done

# Command for all samples (shown but not executed)
: '
for sample in sample*; do
    fastp --in1 ${DATA_DIR}/${sample}_R1.fastq.gz \\
          --in2 ${DATA_DIR}/${sample}_R2.fastq.gz \\
          --out1 Results/fastp/${sample}_R1_cleaned.fastq.gz \\
          --out2 Results/fastp/${sample}_R2_cleaned.fastq.gz \\
          --html Results/fastp/${sample}_report.html \\
          --json Results/fastp/${sample}_report.json \\
          --thread 10
done
'"""
    cells.append(new_code_cell(fastp_code))
    return cells

def create_mapping_cells():
    """Create cells for mapping section"""
    cells = []

    # Add documentation
    doc_md = """## Read Mapping with STAR

### Overview
STAR (Spliced Transcripts Alignment to a Reference) performs RNA-seq alignment:
- Handles splice junctions
- Supports paired-end reads
- Provides mapping statistics

### Important Parameters
- `--runThreadN`: Number of threads
- `--genomeDir`: Path to genome index
- `--readFilesIn`: Input files (R1 R2 for paired-end)
- `--outFileNamePrefix`: Prefix for output files
- `--outSAMtype`: Output format specification

### Documentation
- [STAR Manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)
- [RNA-STAR GitHub](https://github.com/alexdobin/STAR)
"""
    cells.append(new_markdown_cell(doc_md))

    # Add mapping code
    mapping_code = """# Map first two samples with STAR
for sample in sample1 sample2; do
    STAR --runThreadN 10 \\
         --genomeDir reference/star_index \\
         --readFilesIn Results/fastp/${sample}_R1_cleaned.fastq.gz \\
                      Results/fastp/${sample}_R2_cleaned.fastq.gz \\
         --readFilesCommand zcat \\
         --outFileNamePrefix Results/star/${sample}_ \\
         --outSAMtype BAM SortedByCoordinate \\
         --outBAMsortingThreadN 10
done

# Index BAM files
for sample in sample1 sample2; do
    samtools index Results/star/${sample}_Aligned.sortedByCoord.out.bam
done

# Commands for all samples (shown but not executed)
: '
for sample in sample*; do
    STAR --runThreadN 10 \\
         --genomeDir reference/star_index \\
         --readFilesIn Results/fastp/${sample}_R1_cleaned.fastq.gz \\
                      Results/fastp/${sample}_R2_cleaned.fastq.gz \\
         --readFilesCommand zcat \\
         --outFileNamePrefix Results/star/${sample}_ \\
         --outSAMtype BAM SortedByCoordinate \\
         --outBAMsortingThreadN 10

    samtools index Results/star/${sample}_Aligned.sortedByCoord.out.bam
done
'"""
    cells.append(new_code_cell(mapping_code))
    return cells

def create_multiqc_cells():
    """Create cells for MultiQC section"""
    cells = []

    # Add documentation
    doc_md = """## Quality Report Aggregation with MultiQC

### Overview
MultiQC aggregates results from multiple bioinformatics tools:
- FastQC reports
- fastp reports
- STAR alignment logs
- Creates a single comprehensive report

### Important Parameters
- `-o`: Output directory
- `-f`: Force overwrite existing reports

### Documentation
- [MultiQC Documentation](https://multiqc.info/)
- [Available Modules](https://multiqc.info/docs/#available-modules)
"""
    cells.append(new_markdown_cell(doc_md))

    # Add MultiQC code
    multiqc_code = """# Run MultiQC on all results
multiqc -o Results/multiqc Results/"""
    cells.append(new_code_cell(multiqc_code))
    return cells

# Create new notebook
detailed_nb = new_notebook()

# Add title and overview
title_md = """# RNA-seq Data Analysis Pipeline
## Date: Tuesday 26/11/2024

This notebook contains the complete workflow for RNA-seq data analysis, including quality control, preprocessing, and mapping steps. The notebook is designed to be self-explanatory and includes detailed documentation for future reference.

## Overview
This pipeline processes paired-end RNA-seq data from mouse samples. We will:
1. Check the quality of raw reads using FastQC
2. Preprocess reads using fastp
3. Verify the quality of processed reads
4. Download and index the reference genome
5. Map reads to the reference genome using STAR
6. Generate comprehensive quality reports with MultiQC

## Data Organization
- Raw data location: `/srv/data/meg-m2-rnaseq/Data/fastq/raw/`
- Genome annotation: `/srv/data/meg-m2-rnaseq/Genomes/Mmu/GRCm39/extracted/`
- Results will be organized in tool-specific subfolders:
  - FastQC results: `Results/fastqc/`
  - fastp results: `Results/fastp/`
  - STAR mapping results: `Results/star/`
  - MultiQC results: `Results/multiqc/`

## Required Software
- FastQC v0.12.1
- MultiQC v1.13
- fastp v0.23.1
- STAR v2.7.11a
- samtools v1.18

## Important Notes
- This notebook demonstrates the workflow on the first two samples
- Commands for processing all samples are provided as commented code
- Results are organized in tool-specific directories
- The pipeline is designed for mouse RNA-seq data
"""
detailed_nb.cells.append(new_markdown_cell(title_md))

# Add setup code
setup_code = """#!/bin/bash
# Create Results directory structure
mkdir -p Results/{fastqc,fastp,star,multiqc}

# Display the created directory structure
tree Results/"""
detailed_nb.cells.append(new_code_cell(setup_code))

# Add sections
detailed_nb.cells.extend(create_genome_indexing_cells())
detailed_nb.cells.extend(create_quality_check_cells())
detailed_nb.cells.extend(create_preprocessing_cells())
detailed_nb.cells.extend(create_mapping_cells())
detailed_nb.cells.extend(create_multiqc_cells())

# Save the new notebook
with open('detailed_notebook.ipynb', 'w') as f:
    nbformat.write(detailed_nb, f)
