import nbformat
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell, new_raw_cell
import json
import re

def load_notebook(filename):
    with open(filename, 'r', encoding='utf-8') as f:
        return nbformat.read(f, as_version=4)

def create_improved_notebook():
    # Load both notebooks
    pipe06 = load_notebook('Pipe_06-bash_reads-counts-pseudomapping.ipynb')
    ps5 = load_notebook('PS5-ReadCounts-bash-2022-executed.ipynb')

    # Create new notebook
    new_nb = new_notebook()

    # Add title and date
    new_nb.cells.append(new_markdown_cell("""# RNA-seq Data Analysis Pipeline - Part 5: Read Counting
## Tuesday 26/11/2024

This notebook covers read counting and quantification using two complementary approaches:
1. Gene-level quantification with featureCounts
2. Transcript-level quantification with Salmon

<div class="alert alert-info">
<b>Learning Objectives:</b><br>
- Understand different approaches to RNA-seq read counting
- Learn to use featureCounts for gene-level quantification
- Learn to use Salmon for transcript-level quantification
- Interpret read count statistics and quality metrics
</div>
"""))

    # Add system setup with simplified variables
    new_nb.cells.append(new_markdown_cell("""## 1. System Setup
### 1.1 Environment Configuration

We'll set up our working environment with the necessary variables and directory structure."""))

    new_nb.cells.append(new_code_cell("""## Code cell 1 ##
# Base directories
WORK_DIR="/srv/home/${USER}/meg_m2_rnaseq_bash"
DATA_DIR="/srv/data/meg-m2-rnaseq"

# Input/Output directories
RESULTS_DIR="${WORK_DIR}/Results"
BAM_DIR="${RESULTS_DIR}/star"
COUNTS_DIR="${RESULTS_DIR}/counts"
SALMON_DIR="${RESULTS_DIR}/salmon"

# Reference files
GTF_FILE="${DATA_DIR}/Genomes/Mmu/GRCm39/extracted/genome_annotation-M35.gtf"
TRANSCRIPTOME="${DATA_DIR}/Genomes/Mmu/GRCm39/extracted/transcriptome.fa"

# System resources (adjusted for Plasmabio)
CPU_CORES=4
MAX_RAM="6G"

# Tool parameters
STRAND_SPECIFIC=0  # 0=unstranded, 1=stranded, 2=reversely stranded
MIN_MAPPING_QUALITY=10
LIBTYPE="A"  # Automatic detection of library type"""))

    # Add directory creation and validation
    new_nb.cells.append(new_code_cell("""## Code cell 2 ##
# Create output directories
mkdir -p ${COUNTS_DIR} ${SALMON_DIR}

# Validate input files
echo "Checking input files..."
[ -f "${GTF_FILE}" ] && echo "GTF file found" || echo "ERROR: GTF file missing"
[ -d "${BAM_DIR}" ] && echo "BAM directory found" || echo "ERROR: BAM directory missing"
"""))

    # Add featureCounts section
    new_nb.cells.append(new_markdown_cell("""## 2. Gene-Level Read Counting with featureCounts

<div class="alert alert-info">
<b>Tool Information:</b><br>
featureCounts is part of the Subread package and provides fast and accurate read counting for RNA-seq data.
<br><br>
<b>Key Features:</b>
- Supports both single and paired-end reads
- Handles multi-mapping reads
- Provides detailed assignment statistics
- Efficient memory usage
</div>

For more information, visit the [featureCounts documentation](http://subread.sourceforge.net/)."""))

    # Add featureCounts single sample example
    new_nb.cells.append(new_code_cell("""## Code cell 3 ##
# Run featureCounts on a single sample
SAMPLE1=$(ls ${BAM_DIR}/*.bam | head -n 1)
SAMPLE_NAME=$(basename ${SAMPLE1} .bam)

featureCounts \\
    -T ${CPU_CORES} \\
    -s ${STRAND_SPECIFIC} \\
    -Q ${MIN_MAPPING_QUALITY} \\
    -p \\  # paired-end
    -a ${GTF_FILE} \\
    -o ${COUNTS_DIR}/${SAMPLE_NAME}_counts.txt \\
    ${SAMPLE1}

# Display summary
cat ${COUNTS_DIR}/${SAMPLE_NAME}_counts.txt.summary"""))

    # Add section 1.3.2 - Running featureCounts on multiple samples
    new_nb.cells.append(new_markdown_cell("""### 1.3.2 Running featureCounts on Multiple Samples

<div class="alert alert-info">
<b>Important:</b><br>
Processing multiple samples at once is more efficient and ensures consistent parameters across all samples.
</div>

Now we'll process all samples using a loop. This approach allows us to:
1. Process all samples with identical parameters
2. Generate a combined count matrix
3. Save time and reduce potential errors
4. Create consistent output files"""))

    new_nb.cells.append(new_code_cell("""## Code cell 4 ##
# Create a list of all BAM files
BAM_FILES=(${BAM_DIR}/*.bam)
echo "Found ${#BAM_FILES[@]} BAM files"

# Run featureCounts on all samples
featureCounts \\
    -T ${CPU_CORES} \\
    -s ${STRAND_SPECIFIC} \\
    -Q ${MIN_MAPPING_QUALITY} \\
    -p \\
    -a ${GTF_FILE} \\
    -o ${COUNTS_DIR}/all_samples_counts.txt \\
    ${BAM_FILES[@]}

# Create a simplified count matrix
cut -f 1,7- ${COUNTS_DIR}/all_samples_counts.txt | grep -v '^#' > ${COUNTS_DIR}/counts_matrix.txt"""))

    # Add Salmon section
    new_nb.cells.append(new_markdown_cell("""## 2. Pseudo-mapping with Salmon

<div class="alert alert-info">
<b>Tool Information:</b><br>
Salmon performs transcript-level quantification using lightweight algorithms (pseudo-alignment).
<br><br>
<b>Advantages:</b>
- Faster than traditional alignment
- Direct transcript-level quantification
- Bias-aware estimation
- Memory efficient
</div>

For more information, visit the [Salmon documentation](https://salmon.readthedocs.io/)."""))

    # Add Salmon index and quantification
    new_nb.cells.append(new_code_cell("""## Code cell 5 ##
# Index the transcriptome (if not already done)
salmon index \\
    -t ${TRANSCRIPTOME} \\
    -i ${SALMON_DIR}/transcriptome_index \\
    -p ${CPU_CORES}"""))

    new_nb.cells.append(new_code_cell("""## Code cell 6 ##
# Process first two samples
FASTQ_DIR="${DATA_DIR}/fastq/raw"
SAMPLES=($(ls ${FASTQ_DIR}/*_R1.fastq.gz | head -n 2))

for R1 in "${SAMPLES[@]}"; do
    R2=${R1/_R1/_R2}
    SAMPLE=$(basename ${R1} _R1.fastq.gz)

    echo "Processing sample: ${SAMPLE}"

    salmon quant \\
        -i ${SALMON_DIR}/transcriptome_index \\
        -l ${LIBTYPE} \\
        -1 ${R1} \\
        -2 ${R2} \\
        -p ${CPU_CORES} \\
        --validateMappings \\
        -o ${SALMON_DIR}/${SAMPLE}
done"""))

    # Add raw cell for processing all samples
    new_nb.cells.append(new_raw_cell("""## Code cell 7 ##
# Process all samples
for R1 in ${FASTQ_DIR}/*_R1.fastq.gz; do
    R2=${R1/_R1/_R2}
    SAMPLE=$(basename ${R1} _R1.fastq.gz)

    salmon quant \\
        -i ${SALMON_DIR}/transcriptome_index \\
        -l ${LIBTYPE} \\
        -1 ${R1} \\
        -2 ${R2} \\
        -p ${CPU_CORES} \\
        --validateMappings \\
        -o ${SALMON_DIR}/${SAMPLE}
done"""))

    # Add quality control section
    new_nb.cells.append(new_markdown_cell("""## 3. Quality Control

Let's examine the quality metrics for both quantification methods. This helps us ensure the reliability of our results."""))

    new_nb.cells.append(new_code_cell("""## Code cell 8 ##
# Check featureCounts assignment statistics
echo "=== featureCounts Assignment Statistics ==="
cat ${COUNTS_DIR}/all_samples_counts.txt.summary

# Check Salmon mapping rates
echo -e "\\n=== Salmon Mapping Rates ==="
for dir in ${SALMON_DIR}/*/; do
    sample=$(basename ${dir})
    rate=$(grep -A 1 "Mapping Rate" ${dir}/logs/salmon_quant.log | tail -n 1)
    echo "${sample}: ${rate}"
done"""))

    # Add results interpretation section
    new_nb.cells.append(new_markdown_cell("""## 4. Results Interpretation

<div class="alert alert-info">
<b>Key Points:</b><br>
- featureCounts provides gene-level counts
- Salmon provides transcript-level abundance estimates
- Both tools generate quality metrics
- Compare mapping rates between methods
</div>

The output files can be found in:
1. Gene counts: ${COUNTS_DIR}/counts_matrix.txt
2. Transcript quantification: ${SALMON_DIR}/<sample>/quant.sf"""))

    return new_nb

def save_notebook(notebook, filename):
    with open(filename, 'w', encoding='utf-8') as f:
        nbformat.write(notebook, f)

if __name__ == '__main__':
    improved_nb = create_improved_notebook()
    save_notebook(improved_nb, 'PS5-2024-merged.ipynb')
