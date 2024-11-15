import nbformat as nbf

# Create a new notebook
nb = nbf.v4.new_notebook()

# Metadata
nb.metadata = {
    "kernelspec": {
        "display_name": "Bash",
        "language": "bash",
        "name": "bash"
    },
    "language_info": {
        "codemirror_mode": "shell",
        "file_extension": ".sh",
        "mimetype": "text/x-sh",
        "name": "bash"
    }
}

# Title and Overview
title_md = """# Mouse RNA-seq Data Analysis Pipeline
### Master Course - Tuesday 26/11/2024

<div class="alert alert-block alert-info">
<b>Course Overview:</b><br>
This notebook guides you through the first half of the RNA-seq analysis pipeline for mouse paired-end data. You will learn how to:
- Assess raw data quality using FastQC
- Preprocess reads using fastp for adapter removal and quality trimming
- Map cleaned reads to the mouse genome using STAR
- Generate comprehensive quality reports with MultiQC

<b>Resource Allocation:</b><br>
- CPUs: 10 per student
- RAM: 6 GB per student
</div>

<div class="alert alert-block alert-warning">
<b>Important Data Locations:</b><br>
- Raw fastq files: <code>/srv/data/meg-m2-rnaseq/Data/fastq/raw/</code>
- Mouse genome annotation: <code>/srv/data/meg-m2-rnaseq/Genomes/Mmu/GRCm39/extracted/mouse genome_annotation-M35.gtf</code>
- Fastp results: <code>/srv/data/meg-m2-rnaseq/Results/fastp/</code>

<b>Note:</b> For efficiency, we will demonstrate the analysis on the first two samples. Commands for processing all samples are provided as silent cells.
</div>"""

# Create cells list
cells = [
    nbf.v4.new_markdown_cell(title_md),

    # Setup Section
    nbf.v4.new_markdown_cell("""## 1. Environment Setup
<div class="alert alert-block alert-info">
First, we'll check our tool versions and set up our working environment. This ensures reproducibility and proper resource allocation.
</div>"""),

    nbf.v4.new_code_cell("""# Cell 1: Tool Versions
echo "=== Tool Versions ==="
fastqc --version
fastp --version
STAR --version
samtools --version | head -n 2
multiqc --version"""),

    nbf.v4.new_code_cell("""# Cell 2: Directory Setup
# Create results directories
mkdir -p Results/{fastqc,fastp,star,multiqc}
echo "Created analysis directories:"
ls -l Results/

# Set data paths
FASTQ_DIR="/srv/data/meg-m2-rnaseq/Data/fastq/raw"
GENOME_DIR="/srv/data/meg-m2-rnaseq/Genomes/Mmu/GRCm39/extracted"
FASTP_DIR="/srv/data/meg-m2-rnaseq/Results/fastp"

# Set number of threads (80% of available CPUs)
N_THREADS=8  # Based on 10 CPUs per student"""),

    # Quality Check Section
    nbf.v4.new_markdown_cell("""## 2. Raw Data Quality Assessment
<div class="alert alert-block alert-info">
FastQC performs quality control checks on raw sequence data. It provides:
- Basic statistics
- Per base sequence quality
- Per sequence quality scores
- Sequence duplication levels
- Overrepresented sequences
- Adapter content
</div>"""),

    nbf.v4.new_code_cell("""# Cell 3: Initial FastQC Analysis
# Get first two samples
SAMPLES=($(ls ${FASTQ_DIR}/*_R{1,2}.fastq.gz | head -n 4))
echo "Processing samples:"
printf '%s\\n' "${SAMPLES[@]}"

# Run FastQC
for sample in "${SAMPLES[@]}"; do
    echo "Processing $sample..."
    fastqc --outdir Results/fastqc --threads $N_THREADS "$sample"
done"""),

    nbf.v4.new_raw_cell("""# (Cell 4): FastQC Analysis for All Samples
# Process all samples
for sample in ${FASTQ_DIR}/*_R{1,2}.fastq.gz; do
    echo "Processing $sample..."
    fastqc --outdir Results/fastqc --threads $N_THREADS "$sample"
done"""),

    # Preprocessing Section
    nbf.v4.new_markdown_cell("""## 3. Read Preprocessing
<div class="alert alert-block alert-info">
<b>Fastp performs several crucial preprocessing steps:</b>
- Adapter trimming
- Low quality base trimming
- Quality filtering
- Per base quality pruning
- Length filtering

We'll process the first two samples and provide commands for processing all samples.
</div>"""),


    nbf.v4.new_code_cell("""# Cell 5: Fastp Processing (First 2 Samples)
# Process first two paired-end samples
for i in {1..2}; do
    R1=$(ls ${FASTQ_DIR}/*_R1.fastq.gz | head -n $i | tail -n 1)
    R2=${R1/_R1/_R2}
    base=$(basename $R1 _R1.fastq.gz)

    echo "Processing $base..."
    fastp --in1 $R1 \\
          --in2 $R2 \\
          --out1 Results/fastp/${base}_R1.cleaned.fastq.gz \\
          --out2 Results/fastp/${base}_R2.cleaned.fastq.gz \\
          --html Results/fastp/${base}_fastp.html \\
          --json Results/fastp/${base}_fastp.json \\
          --thread $N_THREADS \\
          --detect_adapter_for_pe
done"""),

    nbf.v4.new_raw_cell("""# (Cell 6): Fastp Processing for All Samples
# Process all paired-end samples
for R1 in ${FASTQ_DIR}/*_R1.fastq.gz; do
    R2=${R1/_R1/_R2}
    base=$(basename $R1 _R1.fastq.gz)

    echo "Processing $base..."
    fastp --in1 $R1 \\
          --in2 $R2 \\
          --out1 Results/fastp/${base}_R1.cleaned.fastq.gz \\
          --out2 Results/fastp/${base}_R2.cleaned.fastq.gz \\
          --html Results/fastp/${base}_fastp.html \\
          --json Results/fastp/${base}_fastp.json \\
          --thread $N_THREADS \\
          --detect_adapter_for_pe
done"""),

    # Mapping Section
    nbf.v4.new_markdown_cell("""## 4. Read Mapping
<div class="alert alert-block alert-info">
<b>STAR Mapping Process:</b>
1. First, we need to index the reference genome (provided in a silent cell)
2. Then, we'll map the cleaned reads to the reference
3. Finally, we'll index the resulting BAM files

<b>Note:</b> We'll only process the first two samples. BAM and BAI files for the remaining samples will be provided on the server.
</div>"""),

    nbf.v4.new_raw_cell("""# (Cell 7): STAR Genome Indexing
# Prepare genome index
GENOME_FASTA="${GENOME_DIR}/genome.fa"
GENOME_GTF="${GENOME_DIR}/mouse genome_annotation-M35.gtf"  # Updated path to match requirements
STAR_INDEX="Results/star/genome_index"

mkdir -p $STAR_INDEX

# Generate genome index
STAR --runMode genomeGenerate \\
     --genomeDir $STAR_INDEX \\
     --genomeFastaFiles $GENOME_FASTA \\
     --sjdbGTFfile $GENOME_GTF \\
     --sjdbOverhang 99 \\
     --runThreadN $N_THREADS"""),

    nbf.v4.new_code_cell("""# Cell 8: STAR Mapping (First 2 Samples)
# Map first two samples
for i in {1..2}; do
    R1=$(ls Results/fastp/*_R1.cleaned.fastq.gz | head -n $i | tail -n 1)
    R2=${R1/_R1/_R2}
    base=$(basename $R1 _R1.cleaned.fastq.gz)

    echo "Mapping $base..."
    STAR --genomeDir Results/star/genome_index \\
         --readFilesIn $R1 $R2 \\
         --readFilesCommand zcat \\
         --outFileNamePrefix Results/star/${base}_ \\
         --outSAMtype BAM SortedByCoordinate \\
         --runThreadN $N_THREADS

    # Index BAM file
    samtools index Results/star/${base}_Aligned.sortedByCoord.out.bam
done"""),

    # MultiQC Section
    nbf.v4.new_markdown_cell("""## 5. Quality Reports
<div class="alert alert-block alert-info">
MultiQC aggregates results from:
- FastQC (raw data quality)
- Fastp (preprocessing statistics)
- STAR (mapping metrics)

This provides a comprehensive view of the analysis quality at each step.
</div>"""),

    nbf.v4.new_code_cell("""# Cell 9: MultiQC Report Generation
# Generate comprehensive report for all available results
multiqc --outdir Results/multiqc \\
        --filename "mouse_rnaseq_report" \\
        Results/fastqc \\
        Results/fastp \\
        Results/star""")
]

# Add cells to notebook
nb.cells = cells

# Write the notebook
with open('Mouse_RNAseq_Analysis_2024.ipynb', 'w') as f:
    nbf.write(nb, f)
