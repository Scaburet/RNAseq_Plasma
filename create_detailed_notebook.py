import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
import os, json

def create_genome_indexing_cells():
    """Create cells for genome downloading and indexing section"""
    cells = []

    # Add documentation
    doc_md = """## Reference Genome Download and Indexing

### Overview
The mouse reference genome and its annotation are fundamental for RNA-seq analysis:
- Reference genome: DNA sequence assembly (GRCm39)
- Gene annotation: Gene models and features (GTF format)
- Index: Optimized data structure for rapid read mapping

### Understanding Reference Files
1. **Reference Genome (FASTA)**
   - Contains chromosome sequences
   - Represents mouse genome assembly GRCm39
   - Downloaded from Ensembl
   - Essential for read alignment

2. **Gene Annotation (GTF)**
   - Defines gene structures
   - Contains exon/intron boundaries
   - Identifies gene features
   - Required for splice-aware mapping

3. **STAR Index**
   - Optimized search structure
   - Enables fast read mapping
   - Contains splice junction information
   - Memory-intensive to generate

### Resource Requirements
- Disk Space: ~8GB for genome and index
- Memory: ~32GB recommended for indexing
- CPU: Multi-threaded process (using 10 threads)
- Time: 30-60 minutes depending on resources

### Key Parameters for STAR Indexing
- `--runMode genomeGenerate`: Index generation mode
- `--genomeDir`: Index output directory
- `--genomeFastaFiles`: Reference genome path
- `--sjdbGTFfile`: Gene annotation file
- `--runThreadN`: Number of threads
- `--sjdbOverhang`: Read length - 1 (default: 100)

### Documentation Resources
- [Ensembl Mouse Genome](https://www.ensembl.org/Mus_musculus/Info/Index)
- [STAR Manual - Genome Generation](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf#page=7)
- [GTF Format Specification](https://www.ensembl.org/info/website/upload/gff.html)
- [GRCm39 Assembly Info](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.27)
"""
    cells.append(new_markdown_cell(doc_md))

    # Add genome download code
    download_code = """%%bash
# Download and prepare reference genome
echo "Creating reference directories..."
mkdir -p reference/star_index

echo "Downloading mouse reference genome..."
cd reference
wget https://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz

echo "Creating STAR index..."
STAR --runMode genomeGenerate \\
     --genomeDir star_index \\
     --genomeFastaFiles Mus_musculus.GRCm39.dna.primary_assembly.fa \\
     --sjdbGTFfile "/srv/data/meg-m2-rnaseq/Genomes/Mmu/GRCm39/extracted/mouse genome_annotation-M35.gtf" \\
     --runThreadN 10 \\
     --sjdbOverhang 100

echo "Checking index files..."
ls -lh star_index/
echo "Index generation complete."
cd ..
"""
    cells.append(new_code_cell(download_code, metadata={"tags": ["bash"]}))
    return cells

def create_quality_check_cells():
    """Create cells for quality check section"""
    cells = []

    # Add documentation
    doc_md = """## Quality Control with FastQC

### Overview
FastQC performs quality control checks on raw sequence data. Understanding these quality metrics is crucial for:
- Identifying sequencing problems
- Detecting contamination
- Ensuring data reliability
- Planning preprocessing steps

### Quality Metrics Explained
1. **Basic Statistics**
   - Total sequences
   - Sequence length
   - GC content percentage
   - Important for sample validation

2. **Per Base Sequence Quality**
   - Quality scores across all bases
   - Look for declining quality at read ends
   - Phred scores > 28 are considered good
   - Common drop in quality at read ends

3. **Per Sequence Quality Scores**
   - Overall quality distribution
   - Most sequences should have mean quality > 27
   - Identifies problematic sequences

4. **Sequence Duplication Levels**
   - PCR duplication assessment
   - Expected higher in RNA-seq
   - Helps identify library complexity

5. **Overrepresented Sequences**
   - Adapter contamination check
   - Highly expressed genes
   - rRNA contamination

### Important Parameters
- `-o`: Output directory for reports
- `-t`: Number of threads (we use 10)
- `--noextract`: Keep reports zipped
- `--quiet`: Suppress progress messages

### Documentation Resources
- [FastQC Documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Example Reports](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html)
- [Quality Score Encoding](https://en.wikipedia.org/wiki/FASTQ_format#Quality)
- [RNA-seq QC Best Practices](https://rna-seqblog.com/quality-control-tools/)
"""
    cells.append(new_markdown_cell(doc_md))

    # Add FastQC code
    fastqc_code = """%%bash
# Run FastQC on first two samples
echo "Creating FastQC output directory..."
mkdir -p Results/fastqc

DATA_DIR="/srv/data/meg-m2-rnaseq/Data/fastq/raw"
echo "Running FastQC on first two samples..."
fastqc -o Results/fastqc -t 10 \\
    ${DATA_DIR}/sample1_R1.fastq.gz \\
    ${DATA_DIR}/sample1_R2.fastq.gz \\
    ${DATA_DIR}/sample2_R1.fastq.gz \\
    ${DATA_DIR}/sample2_R2.fastq.gz

echo "FastQC analysis complete. Checking results..."
ls -lh Results/fastqc/

# Command for all samples (shown but not executed)
: '
echo "Running FastQC on all samples..."
fastqc -o Results/fastqc -t 10 ${DATA_DIR}/*.fastq.gz
'"""
    cells.append(new_code_cell(fastqc_code, metadata={"tags": ["bash"]}))
    return cells

def create_preprocessing_cells():
    """Create cells for preprocessing section"""
    cells = []

    # Add documentation
    doc_md = """## Read Preprocessing with fastp

### Overview
fastp performs quality control and preprocessing of sequencing reads. This step is crucial for:
- Removing technical artifacts
- Improving mapping accuracy
- Reducing false positives in downstream analysis
- Generating detailed QC reports

### Understanding the Preprocessing Steps
1. **Adapter Trimming**
   - Removes sequencing adapters automatically
   - Prevents misalignments
   - Improves mapping efficiency
   - Uses overlap analysis for adapter detection

2. **Quality Filtering**
   - Removes low-quality bases from read ends
   - Filters out poor quality reads
   - Uses sliding window approach
   - Configurable quality thresholds

3. **Base Correction (Paired-End Only)**
   - Uses overlapping information
   - Corrects mismatched bases
   - Improves sequence accuracy
   - Only applies to overlapping regions

4. **Length Filtering**
   - Removes too short reads
   - Prevents spurious mappings
   - Configurable minimum length
   - Important for mapping accuracy

### Quality Report Components
1. **Before Processing**
   - Initial quality metrics
   - Base composition
   - Quality score distribution
   - Adapter content

2. **After Processing**
   - Final read statistics
   - Quality improvements
   - Filtering statistics
   - Base quality distribution

### Key Parameters Explained
- `--in1/--in2`: Input paired-end files
- `--out1/--out2`: Output cleaned files
- `--html`: HTML report location
- `--json`: JSON format metrics
- `--thread`: CPU threads (10 recommended)
- `--qualified_quality_phred`: Minimum base quality (default: 15)
- `--length_required`: Minimum read length (default: 36)
- `--cut_front`: Trim low quality bases from read start
- `--cut_tail`: Trim low quality bases from read end
- `--cut_mean_quality`: Average quality threshold

### Documentation Resources
- [fastp GitHub](https://github.com/OpenGene/fastp)
- [Parameters Guide](https://github.com/OpenGene/fastp#all-options)
- [Algorithm Details](https://github.com/OpenGene/fastp#algorithm)
- [Publication](https://doi.org/10.1093/bioinformatics/bty560)
"""
    cells.append(new_markdown_cell(doc_md))

    # Add fastp code
    fastp_code = """%%bash
# Process first two samples with fastp
echo "Creating fastp output directory..."
mkdir -p Results/fastp

DATA_DIR="/srv/data/meg-m2-rnaseq/Data/fastq/raw"
for sample in sample1 sample2; do
    echo "Processing ${sample} with fastp..."

    # Run fastp with detailed parameters
    fastp --in1 ${DATA_DIR}/${sample}_R1.fastq.gz \\
          --in2 ${DATA_DIR}/${sample}_R2.fastq.gz \\
          --out1 Results/fastp/${sample}_R1_cleaned.fastq.gz \\
          --out2 Results/fastp/${sample}_R2_cleaned.fastq.gz \\
          --html Results/fastp/${sample}_report.html \\
          --json Results/fastp/${sample}_report.json \\
          --thread 10 \\
          --qualified_quality_phred 15 \\
          --length_required 36 \\
          --cut_front \\
          --cut_tail \\
          --cut_mean_quality 20

    echo "Checking output files for ${sample}..."
    ls -lh Results/fastp/${sample}*

    echo "Summary of processing for ${sample}:"
    echo "--------------------------------"
    echo "Original files:"
    ls -lh ${DATA_DIR}/${sample}_R*.fastq.gz
    echo "Processed files:"
    ls -lh Results/fastp/${sample}_R*_cleaned.fastq.gz
    echo "Reports:"
    ls -lh Results/fastp/${sample}_report.*
    echo "--------------------------------"
done

# Command for all samples (shown but not executed)
: '
echo "Processing all samples with fastp..."
for sample in sample*; do
    fastp --in1 ${DATA_DIR}/${sample}_R1.fastq.gz \\
          --in2 ${DATA_DIR}/${sample}_R2.fastq.gz \\
          --out1 Results/fastp/${sample}_R1_cleaned.fastq.gz \\
          --out2 Results/fastp/${sample}_R2_cleaned.fastq.gz \\
          --html Results/fastp/${sample}_report.html \\
          --json Results/fastp/${sample}_report.json \\
          --thread 10 \\
          --qualified_quality_phred 15 \\
          --length_required 36 \\
          --cut_front \\
          --cut_tail \\
          --cut_mean_quality 20
done
'"""
    cells.append(new_code_cell(fastp_code, metadata={"tags": ["bash"]}))
    return cells

def create_mapping_cells():
    """Create cells for mapping section"""
    cells = []

    # Add documentation
    doc_md = """## Read Mapping with STAR

### Overview
STAR (Spliced Transcripts Alignment to a Reference) is designed specifically for RNA-seq data:
- Handles splice junctions efficiently
- Supports paired-end reads
- High accuracy and speed
- Generates comprehensive mapping statistics

### Understanding the Mapping Process
1. **Seed Search**
   - Finds Maximum Mappable Prefixes (MMPs)
   - Uses uncompressed suffix arrays
   - Extremely fast and memory-efficient
   - Identifies potential mapping locations

2. **Clustering/Stitching**
   - Clusters and extends seed matches
   - Accounts for sequencing errors
   - Handles splice junctions
   - Identifies chimeric alignments

3. **Scoring and Selection**
   - Evaluates alignment quality
   - Considers mismatches and gaps
   - Handles multi-mappers
   - Assigns mapping quality scores

4. **Output Generation**
   - Creates sorted BAM files
   - Includes mapping statistics
   - Records splice junctions
   - Generates alignment logs

### Key Parameters Explained
- `--runThreadN`: Number of threads (10 recommended)
- `--genomeDir`: Path to STAR index
- `--readFilesIn`: Input FASTQ files (R1 R2 for paired-end)
- `--readFilesCommand`: Command for reading compressed files
- `--outFileNamePrefix`: Prefix for output files
- `--outSAMtype`: Output format (BAM SortedByCoordinate)
- `--outBAMsortingThreadN`: Threads for BAM sorting
- `--limitBAMsortRAM`: Memory limit for sorting

### Understanding Output Files
1. **BAM Files**
   - Contains aligned reads
   - Sorted by coordinate
   - Indexed for quick access
   - Used for downstream analysis

2. **Log Files**
   - Final.out: Summary statistics
   - Log.out: Detailed progress
   - SJ.out.tab: Splice junctions
   - Log.progress.out: Real-time progress

### Documentation Resources
- [STAR Manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)
- [RNA-STAR GitHub](https://github.com/alexdobin/STAR)
- [Mapping Statistics Guide](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/)
- [BAM Format Specification](https://samtools.github.io/hts-specs/SAMv1.pdf)
"""
    cells.append(new_markdown_cell(doc_md))

    # Add mapping code
    mapping_code = """%%bash
# Map first two samples with STAR
echo "Creating STAR output directory..."
mkdir -p Results/star

for sample in sample1 sample2; do
    echo "Mapping ${sample} with STAR..."

    # Create sample-specific output directory
    mkdir -p Results/star/${sample}

    STAR --runThreadN 10 \\
         --genomeDir reference/star_index \\
         --readFilesIn Results/fastp/${sample}_R1_cleaned.fastq.gz \\
                      Results/fastp/${sample}_R2_cleaned.fastq.gz \\
         --readFilesCommand zcat \\
         --outFileNamePrefix Results/star/${sample}/ \\
         --outSAMtype BAM SortedByCoordinate \\
         --outBAMsortingThreadN 10 \\
         --limitBAMsortRAM 31532137230 \\
         --outReadsUnmapped Fastx

    echo "Indexing BAM file for ${sample}..."
    samtools index Results/star/${sample}/Aligned.sortedByCoord.out.bam

    echo "Checking mapping results for ${sample}..."
    ls -lh Results/star/${sample}/
    echo "Mapping statistics:"
    cat Results/star/${sample}/Log.final.out
    echo "--------------------------------"
done

# Commands for all samples (shown but not executed)
: '
echo "Mapping all samples with STAR..."
for sample in sample*; do
    mkdir -p Results/star/${sample}

    STAR --runThreadN 10 \\
         --genomeDir reference/star_index \\
         --readFilesIn Results/fastp/${sample}_R1_cleaned.fastq.gz \\
                      Results/fastp/${sample}_R2_cleaned.fastq.gz \\
         --readFilesCommand zcat \\
         --outFileNamePrefix Results/star/${sample}/ \\
         --outSAMtype BAM SortedByCoordinate \\
         --outBAMsortingThreadN 10 \\
         --limitBAMsortRAM 31532137230 \\
         --outReadsUnmapped Fastx

    samtools index Results/star/${sample}/Aligned.sortedByCoord.out.bam
done
'"""
    cells.append(new_code_cell(mapping_code, metadata={"tags": ["bash"]}))
    return cells

def create_multiqc_cells():
    """Create cells for MultiQC section"""
    cells = []

    # Add documentation
    doc_md = """## Quality Report Aggregation with MultiQC

### Overview
MultiQC aggregates quality control reports from multiple bioinformatics tools into a single comprehensive report:
- Combines results from FastQC, fastp, and STAR
- Creates interactive visualizations
- Enables easy sample comparison
- Identifies systematic issues across samples

### Understanding MultiQC Reports
1. **General Statistics**
   - Sample-level metrics
   - Read counts and mapping rates
   - Quality score distributions
   - Duplication rates

2. **FastQC Integration**
   - Sequence quality histograms
   - GC content distribution
   - Adapter content
   - Overrepresented sequences

3. **fastp Results**
   - Filtering statistics
   - Quality improvements
   - Base quality distributions
   - Read length distribution

4. **STAR Alignment Metrics**
   - Mapping rates
   - Uniquely mapped reads
   - Multi-mapped reads
   - Splice junction statistics

### Report Navigation
1. **Interactive Features**
   - Sortable tables
   - Zoomable plots
   - Exportable data
   - Customizable views

2. **Data Interpretation**
   - Sample comparison
   - Quality threshold indicators
   - Warning flags
   - Failure notifications

### Key Parameters
- `-o`: Output directory location
- `-f`: Force overwrite existing reports
- `-d`: Print directory creation
- `--interactive`: Enable interactive plots
- `--zip-data-dir`: Compress data directory

### Documentation Resources
- [MultiQC Documentation](https://multiqc.info/)
- [Available Modules](https://multiqc.info/docs/#available-modules)
- [Report Examples](https://multiqc.info/examples/rna-seq/multiqc_report.html)
- [Configuration Guide](https://multiqc.info/docs/#configuring-multiqc)
"""
    cells.append(new_markdown_cell(doc_md))

    # Add MultiQC code
    multiqc_code = """%%bash
# Run MultiQC on all results
echo "Creating MultiQC output directory..."
mkdir -p Results/multiqc

echo "Running MultiQC..."
multiqc --force \\
        --interactive \\
        --verbose \\
        --dirs \\
        --dir-pattern '*' \\
        -o Results/multiqc \\
        Results/fastqc Results/fastp Results/star

echo "MultiQC report generation complete."
echo "Report location:"
ls -lh Results/multiqc/multiqc_report.html

echo "Available quality reports:"
echo "------------------------"
echo "FastQC reports:"
ls -l Results/fastqc/*.html
echo "------------------------"
echo "fastp reports:"
ls -l Results/fastp/*.html
echo "------------------------"
echo "STAR logs:"
ls -l Results/star/*/Log.final.out
echo "------------------------"
"""
    cells.append(new_code_cell(multiqc_code, metadata={"tags": ["bash"]}))
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
setup_code = """%%bash
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
    nbf.write(detailed_nb, f)
