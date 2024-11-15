import nbformat
import os
import sys

def read_notebook(file_path):
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            return nbformat.read(f, as_version=4)
    except Exception as e:
        print(f"Error reading notebook {file_path}: {str(e)}")
        sys.exit(1)

def write_notebook(notebook, file_path):
    try:
        with open(file_path, 'w', encoding='utf-8') as f:
            nbformat.write(notebook, f)
    except Exception as e:
        print(f"Error writing notebook {file_path}: {str(e)}")
        sys.exit(1)

def create_markdown_cell(content):
    return nbformat.v4.new_markdown_cell(content)

def enhance_notebook():
    # Read both notebooks
    pipe05_path = "/home/ubuntu/workspace/rnaseq/selected_notebooks/Pipe_05-bash_mapping-quality.ipynb"
    merged_path = "/home/ubuntu/workspace/rnaseq/CEA_rnaseq/temp_repo/PS4-2024-merged.ipynb"

    pipe05 = read_notebook(pipe05_path)
    merged = read_notebook(merged_path)

    # Enhanced tool documentation with more beginner-friendly explanations
    tool_docs = """
# RNA-seq Data Quality Assessment Tools

## Understanding the Tools We Use

### 1. FastQC
FastQC helps us check the quality of our sequencing data. Think of it as a quality inspector for our raw data.
- **What it does**: Examines your sequencing files and creates easy-to-read reports
- **Why we use it**: To identify potential problems in our data before we start analysis
- **Key things it checks**:
  - Quality scores for each base
  - GC content (should match what we expect for our organism)
  - Overrepresented sequences (might be contamination or adapters)

[📚 FastQC Documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

### 2. MultiQC
MultiQC takes all our different quality reports and combines them into one easy-to-read report.
- **What it does**: Combines reports from different tools into one
- **Why we use it**: Makes it easier to compare results across samples
- **Input**: Results from FastQC, STAR, QUALIMAP, etc.

[📚 MultiQC Documentation](https://multiqc.info/)

### 3. QUALIMAP
QUALIMAP helps us check how well our reads mapped to the reference genome.
- **What it does**: Analyzes alignment data quality
- **Why we use it**: To ensure our mapping step worked correctly
- **Key metrics it provides**:
  - How many reads mapped to the genome
  - How evenly the reads are distributed
  - Quality of the mapping

[📚 QUALIMAP Documentation](http://qualimap.conesalab.org/)

## Important Parameters Explained

### STAR Alignment Parameters
Let's understand each parameter we use:

```bash
STAR --runMode genomeGenerate \\
     --genomeDir /path/to/genome \\
     --genomeFastaFiles genome.fa \\
     --sjdbGTFfile annotation.gtf \\
     --runThreadN 8 \\
     --sjdbOverhang 99
```

- `--runMode genomeGenerate`: Tells STAR to create an index of the genome
  - Think of this like creating an index in a book to find things quickly
- `--genomeDir`: Where to save the genome index
- `--genomeFastaFiles`: The reference genome file
- `--sjdbGTFfile`: Gene annotation file that tells us where genes are located
- `--runThreadN`: Number of CPU cores to use (speeds up the process)
- `--sjdbOverhang`: Related to read length, usually (read length - 1)

### Samtools Parameters
Samtools helps us work with alignment files. Here are the key commands:

```bash
# Sort BAM file
samtools sort input.bam -o sorted.bam

# Create index
samtools index sorted.bam
```

- `sort`: Organizes the alignments by position in the genome
- `index`: Creates an index for quick access to specific regions
- `-o`: Specifies the output file name
"""

    # Add detailed workflow explanation
    workflow_docs = """
# RNA-seq Analysis Workflow: Step by Step Guide

## Overview
RNA sequencing (RNA-seq) helps us understand which genes are active in our samples. Let's break down each step of our analysis:

## 1. Quality Control of Raw Data
First, we need to check if our sequencing data is good quality:
```bash
# Create directory for FastQC results
mkdir -p Results/fastqc

# Run FastQC on our first two samples
fastqc -o Results/fastqc sample1_R1.fastq.gz sample1_R2.fastq.gz
fastqc -o Results/fastqc sample2_R1.fastq.gz sample2_R2.fastq.gz
```
- `-o`: Specifies output directory
- Why two files per sample? These are paired-end reads (forward and reverse)

## 2. Read Preprocessing
Clean up our raw reads to remove low-quality sequences:
```bash
# Create directory for cleaned reads
mkdir -p Results/fastp

# Process first sample
fastp -i sample1_R1.fastq.gz -I sample1_R2.fastq.gz \\
      -o Results/fastp/sample1_R1_clean.fastq.gz \\
      -O Results/fastp/sample1_R2_clean.fastq.gz \\
      --html Results/fastp/sample1_report.html
```
- `-i/-I`: Input files (forward/reverse reads)
- `-o/-O`: Output files (cleaned reads)
- `--html`: Report file

## 3. Read Mapping with STAR
Map our cleaned reads to the reference genome:
```bash
# Create output directory
mkdir -p Results/star

# Map reads for first sample
STAR --genomeDir reference/star_index \\
     --readFilesIn Results/fastp/sample1_R1_clean.fastq.gz \\
                   Results/fastp/sample1_R2_clean.fastq.gz \\
     --outFileNamePrefix Results/star/sample1_ \\
     --runThreadN 8 \\
     --outSAMtype BAM SortedByCoordinate
```

## 4. Quality Assessment of Mapped Reads
Check how well our reads mapped:
```bash
# Create directory for QUALIMAP results
mkdir -p Results/qualimap/sample1

# Run QUALIMAP
qualimap rnaseq -bam Results/star/sample1_Aligned.sortedByCoord.out.bam \\
                -gtf reference/annotation.gtf \\
                -outdir Results/qualimap/sample1 \\
                --java-mem-size=6G
```

## Directory Structure
```
Results/
├── fastqc/          # FastQC quality reports
├── fastp/           # Cleaned reads and reports
├── star/            # Mapped reads (BAM files)
├── qualimap/        # Mapping quality reports
│   └── sample1/
└── multiqc/         # Combined quality report
```

## Important Notes for Students
- Always check tool versions before starting
- Keep your files organized in appropriate directories
- Monitor your computer's resources (CPU, memory)
- Save your commands in a script for reproducibility
- Take notes on parameter choices and their effects
"""

    # Add progressive loop building explanation
    loop_docs = """
# Understanding and Building Loops for Multiple Samples

## Why Use Loops?
When we have multiple samples, typing the same commands repeatedly is:
1. Time-consuming
2. Error-prone
3. Not efficient

Let's learn how to automate this process step by step!

## Building a Loop for QUALIMAP

### Step 1: List Our Samples
First, let's see how we would run QUALIMAP on one sample:
```bash
qualimap rnaseq -bam Results/star/sample1_Aligned.sortedByCoord.out.bam \\
                -gtf reference/annotation.gtf \\
                -outdir Results/qualimap/sample1 \\
                --java-mem-size=6G
```

### Step 2: Identify the Parts That Change
For each sample, we need to change:
1. Input BAM file name
2. Output directory name

### Step 3: Create a Simple Loop
Let's start with a basic loop for two samples:
```bash
for sample in sample1 sample2; do
    echo "Processing $sample..."
    qualimap rnaseq -bam Results/star/${sample}_Aligned.sortedByCoord.out.bam \\
                    -gtf reference/annotation.gtf \\
                    -outdir Results/qualimap/${sample} \\
                    --java-mem-size=6G
done
```

### Step 4: Make It More Flexible
Now let's modify our loop to handle any number of samples:
```bash
# Get list of all BAM files
BAM_FILES=$(ls Results/star/*_Aligned.sortedByCoord.out.bam)

# Process each BAM file
for bam in $BAM_FILES; do
    # Extract sample name from BAM file
    sample=$(basename $bam _Aligned.sortedByCoord.out.bam)

    echo "Processing sample: $sample"

    # Create output directory
    mkdir -p Results/qualimap/${sample}

    # Run QUALIMAP
    qualimap rnaseq -bam $bam \\
                    -gtf reference/annotation.gtf \\
                    -outdir Results/qualimap/${sample} \\
                    --java-mem-size=6G
done
```

### Understanding Each Part
- `$(ls Results/star/*_Aligned.sortedByCoord.out.bam)`: Lists all BAM files
- `basename`: Extracts the sample name from the file path
- `mkdir -p`: Creates the output directory if it doesn't exist

### Tips for Success
1. Always test your loop with `echo` first
2. Start with a few samples before running on all
3. Monitor the progress with echo statements
4. Keep your file naming consistent
"""

    # Insert new cells at appropriate positions
    new_cells = []
    for cell in merged.cells:
        new_cells.append(cell)
        if "## Tools Overview" in str(cell.source):
            new_cells.append(create_markdown_cell(tool_docs))
        elif "## Workflow" in str(cell.source):
            new_cells.append(create_markdown_cell(workflow_docs))
        elif "## Loop" in str(cell.source):
            new_cells.append(create_markdown_cell(loop_docs))

    merged.cells = new_cells

    # Write enhanced notebook
    write_notebook(merged, merged_path)
    print("Successfully enhanced notebook with detailed explanations")

if __name__ == "__main__":
    enhance_notebook()
