# Summary of Pipe_06-bash_reads-counts.ipynb

## Purpose
This notebook performs RNA sequencing read counting using two complementary approaches:
1. Gene-level quantification with featureCounts
2. Transcript-level quantification with Salmon

## Main Components

### 1. Initial Setup
- Resource configuration (CPU and RAM allocation)
- Directory structure setup
- Input validation

### 2. featureCounts Analysis
- BAM file processing
- Gene-level read counting
- Strand-specific counting options
- Output format generation

### 3. Salmon Analysis
- Transcript-level quantification
- Pseudo-alignment approach
- Expression level estimation

### 4. Quality Control
- Count statistics generation
- Read assignment categorization
- Summary report creation

## Educational Value
- Understanding read counting methods
- Learning different quantification approaches
- Hands-on experience with industry tools
- Data interpretation skills

## Technical Requirements
- Input: Sorted BAM files
- Tools: featureCounts, Salmon
- Resources: 4 CPU cores, 20GB RAM
- Reference: GTF annotation file

## Output Files
- Gene-level count matrices
- Transcript abundance estimates
- Quality control reports
- Assignment statistics
