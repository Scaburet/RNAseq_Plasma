import nbformat
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell, new_raw_cell
import datetime

def create_merged_notebook():
    """Create a merged notebook with intermediate verbosity level."""
    nb = new_notebook()

    # Header and Introduction
    nb.cells.extend([
        new_markdown_cell("""# RNA-seq Analysis Module - Mapping Quality Assessment
## Practical Session 4 (Tuesday 26/11/2024)
---
This notebook guides you through the quality assessment of RNA-seq mapping outputs, focusing on:
- BAM file quality metrics using samtools
- Detailed mapping quality analysis with QUALIMAP
- Multi-sample quality report generation with MultiQC

<div class="alert alert-block alert-info">
<b>Prerequisites:</b><br>
- Completed mapping of RNA-seq reads to reference genome
- Sorted BAM files (.bam) and their indexes (.bai)
- Reference genome annotation file (GTF format)
</div>

<div class="alert alert-block alert-warning">
<b>Resource Requirements:</b><br>
This analysis requires:
- Approximately 30-45 minutes to complete
- 10 CPU cores maximum
- 6GB RAM per student
- ~20GB disk space for results
</div>
"""),

        # Resource Setup
        new_markdown_cell("""## 1. Environment Setup
### 1.1 Resource Configuration

Before starting the analysis, we need to configure our environment and set resource limits. This ensures efficient processing while staying within the available computational resources.

Key parameters:
- NCPUS: Maximum number of CPU cores to use
- MAXRAM: Maximum RAM allocation per process
"""),

        new_code_cell("""## Code Cell n°1 ##
# Maximum resources available
NCPUS=10
MAXRAM="6G"  # Maximum RAM per student
echo "Using maximum ${NCPUS} CPUs and ${MAXRAM} RAM"
"""),

        # Directory Setup
        new_markdown_cell("""### 1.2 Directory Structure

The analysis requires a specific directory structure to organize input data and results. We'll create separate directories for each tool's output to maintain clarity and organization.

<div class="alert alert-block alert-info">
<b>Directory Organization:</b><br>
- Results/samtools/: Basic BAM statistics and metrics
- Results/qualimap/: Detailed mapping quality analysis
- Results/multiqc/: Combined quality reports
</div>
"""),

        new_code_cell("""## Code Cell n°2 ##
# Create Results directory structure
mkdir -p Results/{samtools,qualimap,multiqc}

# Define data paths
DATA_DIR="/srv/data/meg-m2-rnaseq/Data"
FASTQ_DIR="${DATA_DIR}/fastq/raw"
GENOME_DIR="/srv/data/meg-m2-rnaseq/Genomes/Mmu/GRCm39/extracted"
GTF_FILE="${GENOME_DIR}/genome_annotation-M35.gtf"

echo "Results directories created"
ls -l Results/
"""),

        # Reference Genome Section
        new_markdown_cell("""## 2. Reference Genome Preparation
### 2.1 Genome Files

For quality assessment, we need the reference genome and its annotation. These files are used by QUALIMAP to analyze gene coverage and other metrics.

<div class="alert alert-block alert-info">
<b>Required Files:</b><br>
- Genome sequence (FASTA format): Contains the reference genome sequence
- Gene annotation (GTF format): Contains gene structure information

<b>Documentation:</b><br>
- [Ensembl Mouse Genome](http://www.ensembl.org/Mus_musculus/Info/Index)
- [GRCm39 Assembly](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.27)
</div>

<div class="alert alert-block alert-warning">
<b>Important Notes:</b><br>
- Mouse genome version: GRCm39 (Ensembl release 109)
- Files are large (>1GB), download may take time
- Genome and annotation versions must match
</div>

The following code shows how to download and verify the reference files:"""),

        new_raw_cell("""## Code Cell n°(3) ##
# Download reference files if needed
mkdir -p ${GENOME_DIR}
wget -P ${GENOME_DIR} ftp://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
wget -P ${GENOME_DIR} ftp://ftp.ensembl.org/pub/release-109/gtf/mus_musculus/Mus_musculus.GRCm39.109.gtf.gz

# Extract files
gunzip ${GENOME_DIR}/*.gz

# Verify downloads
ls -lh ${GENOME_DIR}/*.{fa,gtf}
"""),

        new_code_cell("""## Code Cell n°3 ##
# Verify reference files
echo "=== Checking genome file ==="
head -n 2 ${GENOME_DIR}/Mus_musculus.GRCm39.dna.primary_assembly.fa
echo -e "\\nChromosome count:"
grep "^>" ${GENOME_DIR}/Mus_musculus.GRCm39.dna.primary_assembly.fa | wc -l

echo -e "\\n=== Checking annotation file ==="
head -n 2 ${GENOME_DIR}/genome_annotation-M35.gtf
echo -e "\\nProtein-coding genes:"
grep -c "gene_biotype \\"protein_coding\\"" ${GENOME_DIR}/genome_annotation-M35.gtf
"""),

        # BAM File Analysis
        new_markdown_cell("""## 3. BAM File Quality Assessment
### 3.1 Initial BAM File Inspection

BAM files contain aligned sequencing reads in a binary format. We'll use samtools to examine these files and generate basic statistics.

<div class="alert alert-block alert-info">
<b>Samtools Functions:</b><br>
- view: Convert between SAM/BAM formats and inspect alignments
- sort: Sort alignments by coordinates or read names
- index: Create index for fast random access
- stats: Generate comprehensive mapping statistics

<b>Key Statistics:</b><br>
- Total reads: Number of sequenced fragments
- Mapped reads: Successfully aligned to reference
- Properly paired: Both mates aligned correctly
- Insert size: Distance between paired reads
</div>

<div class="alert alert-block alert-warning">
<b>Important Parameters:</b><br>
- -F/-f: Include/exclude reads based on flags
- -q: Minimum mapping quality score
- -@ threads: Number of threads to use
</div>

<div class="alert alert-block alert-info">
<b>Documentation:</b><br>
- [Samtools Documentation](http://www.htslib.org/doc/samtools.html)
- [SAM Format Specification](https://samtools.github.io/hts-specs/SAMv1.pdf)
</div>

First, let's examine our BAM files and generate basic statistics:"""),


        new_code_cell("""## Code Cell n°4 ##
# List available BAM files
echo "Available BAM files:"
ls -lh Results/star/*.bam

# Generate basic statistics for first two samples
for bamfile in $(ls Results/star/*.bam | head -n 2); do
    echo "Processing ${bamfile}..."
    sample=$(basename ${bamfile} _Aligned.sortedByCoord.out.bam)

    # Create output directory
    mkdir -p Results/samtools/${sample}

    # Basic alignment statistics
    echo "=== Basic Statistics for ${sample} ==="
    samtools flagstat ${bamfile} > Results/samtools/${sample}/flagstat.txt
    cat Results/samtools/${sample}/flagstat.txt

    # Detailed statistics
    echo -e "\\n=== Detailed Mapping Statistics ==="
    samtools stats ${bamfile} > Results/samtools/${sample}/stats.txt

    # Show summary metrics
    echo "Summary statistics:"
    grep ^SN Results/samtools/${sample}/stats.txt | cut -f 2-

    # Generate insert size metrics
    echo -e "\\n=== Insert Size Distribution ==="
    samtools stats ${bamfile} | grep ^IS | cut -f 2- > Results/samtools/${sample}/insert_size.txt
    head -n 20 Results/samtools/${sample}/insert_size.txt

    # Check for unmapped reads
    echo -e "\\n=== Unmapped Reads Analysis ==="
    samtools view -c -f 4 ${bamfile}

    # Check for quality distribution
    echo -e "\\n=== Mapping Quality Distribution ==="
    samtools view ${bamfile} | cut -f5 | sort -n | uniq -c | head -n 10
done
"""),

        new_markdown_cell("""### 3.2 Understanding BAM Statistics

The statistics generated by samtools provide important insights into the quality of your RNA-seq data:

1. **Flagstat Output:**
   - Total reads: All sequences in the BAM file
   - Mapped reads: Successfully aligned to the reference
   - Properly paired: Both mates aligned as expected
   - Secondary alignments: Alternative alignments for the same read
   - Duplicates: PCR or optical duplicates
   - Supplementary alignments: Split alignments

2. **Insert Size Distribution:**
   - Mean and standard deviation
   - Peak insert size frequency
   - Distribution shape (normal vs. skewed)
   - Important for paired-end data quality assessment

3. **Mapping Quality:**
   - MAPQ scores range from 0 to 60
   - Higher scores indicate more unique alignments
   - Score of 0 indicates multiple mapping locations

<div class="alert alert-block alert-success">
<b>Quality Indicators:</b><br>
✓ High mapping rate (>80%)
✓ High proper pair rate (>75%)
✓ Consistent insert size distribution
✓ High proportion of unique alignments
✗ High duplicate rate (>20%)
✗ Many unmapped reads
✗ High proportion of MAPQ=0 reads
</div>

<div class="alert alert-block alert-info">
<b>Troubleshooting:</b><br>
- Low mapping rate: Check for contamination or wrong reference genome
- High duplicate rate: Review PCR amplification steps
- Inconsistent insert size: Review library preparation protocol
- Many MAPQ=0 reads: Consider using more stringent mapping parameters
</div>"""),

        new_raw_cell("""## Code Cell n°(5) ##
# Process all BAM files
for bamfile in Results/star/*.bam; do
    echo "Processing ${bamfile}..."
    sample=$(basename ${bamfile} _Aligned.sortedByCoord.out.bam)

    # Create output directory
    mkdir -p Results/samtools/${sample}

    # Generate all statistics
    samtools flagstat ${bamfile} > Results/samtools/${sample}/flagstat.txt
    samtools stats ${bamfile} > Results/samtools/${sample}/stats.txt

    # Generate insert size metrics
    samtools stats ${bamfile} | grep ^IS | cut -f 2- > Results/samtools/${sample}/insert_size.txt

    # Generate mapping quality distribution
    samtools view ${bamfile} | cut -f5 | sort -n | uniq -c > Results/samtools/${sample}/mapq_dist.txt
done
"""),

        # QUALIMAP Analysis
        new_markdown_cell("""## 4. QUALIMAP Analysis
### 4.1 Tool Introduction

QUALIMAP is a platform-independent application that provides comprehensive quality control analysis of alignment sequencing data.

Key features:
- Multi-sample BAM QC analysis
- RNA-seq specific quality metrics
- Detailed coverage analysis
- HTML report generation

<div class="alert alert-block alert-warning">
<b>Important Parameters:</b><br>
- -bam: Input BAM file
- -c: Collect overlap rate of reads
- -gd: Genome distribution computation
- -outdir: Output directory
- -nt: Number of threads
- --java-mem-size: Maximum memory allocation
</div>

<div class="alert alert-block alert-info">
<b>Documentation:</b><br>
For more information, visit:
- [QUALIMAP Documentation](http://qualimap.conesalab.org/)
- [QUALIMAP RNA-seq QC](http://qualimap.conesalab.org/doc_html/analysis.html#rna-seq-qc)
</div>"""),

        new_code_cell("""## Code Cell n°6 ##
# Create output directory for QUALIMAP results
mkdir -p Results/qualimap/{bamqc,rnaseq}

# Run QUALIMAP bamqc on first two samples
for bamfile in $(ls Results/star/*.bam | head -n 2); do
    sample=$(basename ${bamfile} _Aligned.sortedByCoord.out.bam)
    echo "Processing ${sample}..."

    qualimap bamqc \
        -bam ${bamfile} \
        -c \
        -gd MOUSE \
        -outdir Results/qualimap/bamqc/${sample} \
        -nt ${NCPUS} \
        --java-mem-size=${MAXRAM}
done

# Run RNA-seq specific analysis
for bamfile in $(ls Results/star/*.bam | head -n 2); do
    sample=$(basename ${bamfile} _Aligned.sortedByCoord.out.bam)
    echo "Processing RNA-seq QC for ${sample}..."

    qualimap rnaseq \
        -bam ${bamfile} \
        -gtf ${GTF_FILE} \
        --paired \
        -outdir Results/qualimap/rnaseq/${sample} \
        -pe \
        --java-mem-size=${MAXRAM}
done
"""),

        new_raw_cell("""## Code Cell n°(7) ##
# Process all samples with QUALIMAP
for bamfile in Results/star/*.bam; do
    sample=$(basename ${bamfile} _Aligned.sortedByCoord.out.bam)

    # Run bamqc
    qualimap bamqc \
        -bam ${bamfile} \
        -c \
        -gd MOUSE \
        -outdir Results/qualimap/bamqc/${sample} \
        -nt ${NCPUS} \
        --java-mem-size=${MAXRAM}

    # Run RNA-seq QC
    qualimap rnaseq \
        -bam ${bamfile} \
        -gtf ${GTF_FILE} \
        --paired \
        -outdir Results/qualimap/rnaseq/${sample} \
        -pe \
        --java-mem-size=${MAXRAM}
done
"""),

        # MultiQC Analysis
        new_markdown_cell("""## 5. MultiQC Report Generation
### 5.1 Tool Overview

MultiQC aggregates results from bioinformatics analyses across many samples into a single report. It automatically scans given directories for recognized log files and compiles a HTML report with interactive plots.

<div class="alert alert-block alert-info">
<b>Key Features:</b>
- Combines QC reports from multiple tools
- Interactive plots and tables
- Easy to interpret summary statistics
- Exportable plots and data

<b>Documentation:</b><br>
For more information, visit:
- [MultiQC Documentation](https://multiqc.info/)
- [Available MultiQC Modules](https://multiqc.info/docs/#multiqc-modules)
</div>
"""),

        new_code_cell("""## Code Cell n°8 ##
# Set MultiQC parameters
MULTIQC_TITLE="Mouse RNA-seq Quality Report"
MULTIQC_COMMENT="Quality control analysis of mouse RNA-seq data"
MULTIQC_FILENAME="mouse_rnaseq_multiqc"

# Run MultiQC
multiqc \
    Results/{samtools,qualimap}/ \
    --title "${MULTIQC_TITLE}" \
    --comment "${MULTIQC_COMMENT}" \
    --filename ${MULTIQC_FILENAME} \
    --outdir Results/multiqc/
"""),

        # Results Interpretation
        new_markdown_cell("""## 6. Results Interpretation

After running the quality control pipeline, you should examine several key metrics:

1. **Basic Mapping Statistics** (samtools):
   - Total reads mapped
   - Properly paired reads
   - Insert size distribution

2. **QUALIMAP Metrics**:
   - Coverage distribution
   - Gene body coverage
   - Read genomic origin
   - Transcript coverage uniformity

3. **MultiQC Summary**:
   - Sample comparison
   - Quality trends
   - Potential batch effects

<div class="alert alert-block alert-success">
<b>Questions to Consider:</b><br>
1. Are the mapping rates consistent across samples?
2. How uniform is the coverage across gene bodies?
3. Are there any concerning quality metrics in specific samples?
4. How do the samples compare in terms of sequencing depth?
</div>
"""),

        # Conclusion
        new_markdown_cell("""## 7. Analysis Complete
---
<div class="alert alert-block alert-success">
<b>Completion:</b><br>
You have successfully:
- Examined BAM file quality using samtools
- Generated detailed quality metrics with QUALIMAP
- Created a comprehensive MultiQC report

The results are stored in:
- Results/samtools/: Basic BAM statistics
- Results/qualimap/: Detailed mapping quality analysis
- Results/multiqc/: Combined quality report

<b>Next Steps:</b><br>
- Review the MultiQC report for overall quality assessment
- Investigate any samples with unusual metrics
- Document any quality concerns for downstream analysis
</div>
""")
    ])

    return nb

# Create and save the notebook
nb = create_merged_notebook()
output_file = "PS4-2024-merged.ipynb"
with open(output_file, 'w') as f:
    nbformat.write(nb, f)

print(f"Created notebook: {output_file}")
