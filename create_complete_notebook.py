import nbformat
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell, new_raw_cell

def create_notebook():
    nb = new_notebook()

    # Title and Introduction
    nb.cells.append(new_markdown_cell("""# RNA-seq Data Analysis Pipeline - PS4 2024
## Tuesday 26/11/2024

This notebook covers the analysis of RNA-seq data from mouse samples, focusing on quality assessment and processing of mapped reads.

<div class="alert alert-info">
<b>Environment Requirements:</b><br>
- 10 CPUs per student<br>
- 6 GB RAM per student<br>
- Tools: samtools, QUALIMAP, MultiQC
</div>
"""))

    # Required Files Section
    nb.cells.append(new_markdown_cell("""## Required Files and Data

<div class="alert alert-warning">
<b>Required Input Files:</b>
1. Paired-end RNA-seq data (fastq.gz format)
   - Location: `/srv/data/meg-m2-rnaseq/Data/fastq/raw/`
   - Format: `*_1.fastq.gz` and `*_2.fastq.gz` for paired-end reads

2. Reference Genome Files:
   - Mouse genome (GRCm39)
   - GTF annotation: `/srv/data/meg-m2-rnaseq/Genomes/Mmu/GRCm39/extracted/genome_annotation-M35.gtf`

3. Previously Generated Files:
   - BAM files from STAR mapping
   - fastp results in `/srv/data/meg-m2-rnaseq/Results/fastp/`
</div>

<div class="alert alert-info">
<b>Output Structure:</b>
Results will be organized in tool-specific directories:
- `Results/samtools/` - Alignment statistics
- `Results/qualimap/` - RNA-seq QC metrics
- `Results/multiqc/` - Combined quality reports
</div>
"""))

    # Reference Genome Section
    nb.cells.append(new_markdown_cell("""## Reference Genome and Annotation

This section covers the downloading and indexing of reference genome files. These steps are typically performed during environment setup.

<div class="alert alert-info">
<b>Reference Files:</b>
- Genome: Mus musculus GRCm39
- Source: Ensembl Release 109
- Annotation: genome_annotation-M35.gtf
</div>

<div class="alert alert-warning">
<b>Note:</b> The following commands are provided for reference and documentation. They are kept in raw format as they are typically executed during environment setup.
</div>
"""))

    # Genome Download and Index Commands
    nb.cells.append(new_raw_cell("""# Create genome directory
mkdir -p reference/mouse/GRCm39

# Download mouse reference genome
wget -P reference/mouse/GRCm39 https://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz

# Extract genome
gunzip reference/mouse/GRCm39/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz

# Index the genome
cd reference/mouse/GRCm39
samtools faidx Mus_musculus.GRCm39.dna.primary_assembly.fa"""))

    # Quality Assessment Introduction
    nb.cells.append(new_markdown_cell("""## Quality Assessment of Mapped Reads

We'll use multiple tools to assess the quality of our mapped reads:

1. **samtools**: Basic alignment statistics
   - Read counts
   - Mapping quality
   - Insert size distribution

2. **QUALIMAP**: Detailed RNA-seq metrics
   - Gene coverage
   - Read distribution
   - Transcript coverage

3. **MultiQC**: Combined report generation
   - Aggregates results from all tools
   - Creates interactive visualizations

<div class="alert alert-info">
<b>Documentation Links:</b>
- [Samtools Manual](http://www.htslib.org/doc/samtools.html)
- [QUALIMAP Documentation](http://qualimap.conesalab.org/doc_html/index.html)
- [MultiQC Documentation](https://multiqc.info/)
</div>
"""))

    # Directory Setup
    nb.cells.append(new_markdown_cell("""### Setup Results Directories

First, we'll create the necessary directories for organizing our results.
"""))

    nb.cells.append(new_code_cell("""# Create results directories
mkdir -p Results/samtools Results/qualimap Results/multiqc"""))

    # Process First Two Samples
    nb.cells.append(new_markdown_cell("""### Process First Two Samples

We'll demonstrate the analysis pipeline on the first two samples. The same process will be applied to all samples later.

<div class="alert alert-info">
<b>Important Parameters:</b>
- samtools stats: Generates comprehensive alignment statistics
- QUALIMAP rnaseq:
  - java-mem-size=6G: Memory allocation
  - GTF file: Required for gene-based analysis
</div>
"""))

    nb.cells.append(new_code_cell("""# Get first two samples
samples=$(ls /srv/data/meg-m2-rnaseq/Data/fastq/raw/*_1.fastq.gz | head -n 2)

# Process each sample
for sample in $samples; do
    # Extract base name
    base=$(basename $sample _1.fastq.gz)
    echo "Processing sample: $base"

    # Samtools stats
    echo "Running samtools stats..."
    samtools stats ${base}.bam > Results/samtools/${base}_stats.txt

    # QUALIMAP analysis
    echo "Running QUALIMAP..."
    qualimap rnaseq \\
        -bam ${base}.bam \\
        -gtf /srv/data/meg-m2-rnaseq/Genomes/Mmu/GRCm39/extracted/genome_annotation-M35.gtf \\
        --java-mem-size=6G \\
        -outdir Results/qualimap/${base}
done"""))

    # Commands for All Samples
    nb.cells.append(new_markdown_cell("""### Process All Samples

The following commands (in raw format) show how to process all samples in the dataset.
"""))

    nb.cells.append(new_raw_cell("""# Process all samples
for sample in /srv/data/meg-m2-rnaseq/Data/fastq/raw/*_1.fastq.gz; do
    base=$(basename $sample _1.fastq.gz)
    echo "Processing sample: $base"

    # Samtools stats
    samtools stats ${base}.bam > Results/samtools/${base}_stats.txt

    # QUALIMAP analysis
    qualimap rnaseq \\
        -bam ${base}.bam \\
        -gtf /srv/data/meg-m2-rnaseq/Genomes/Mmu/GRCm39/extracted/genome_annotation-M35.gtf \\
        --java-mem-size=6G \\
        -outdir Results/qualimap/${base}
done"""))

    # MultiQC Section
    nb.cells.append(new_markdown_cell("""## MultiQC Report Generation

MultiQC combines the quality reports from all samples into a comprehensive report.

<div class="alert alert-info">
<b>MultiQC Features:</b>
- Combines reports from multiple tools
- Creates interactive visualizations
- Enables easy sample comparison
- Generates HTML report

<b>Input Sources:</b>
- samtools stats results
- QUALIMAP RNA-seq reports
</div>
"""))

    nb.cells.append(new_code_cell("""# Generate MultiQC report
echo "Generating MultiQC report..."
multiqc \\
    Results/samtools/ \\
    Results/qualimap/ \\
    -o Results/multiqc/ \\
    --title "Mouse RNA-seq Quality Report" \\
    --comment "Generated for PS4-2024 course"
"""))

    # Save notebook
    with open('PS4-2024-merged.ipynb', 'w', encoding='utf-8') as f:
        nbformat.write(nb, f)

if __name__ == '__main__':
    create_notebook()
