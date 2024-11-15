import nbformat as nbf

# Create a new notebook
nb = nbf.v4.new_notebook()

# Title and overview
title_md = """# RNA-seq Data Analysis Pipeline
### Master Course - Tuesday 26/11/2024

<div class="alert alert-info">
<b>Course Overview:</b><br>
This notebook covers the essential steps of RNA-seq data analysis including:
- Quality control of raw sequencing data
- Read preprocessing and quality verification
- Read mapping to reference genome
</div>

<div class="alert alert-warning">
<b>Resource Requirements:</b><br>
- CPUs: 10
- RAM: 6GB
</div>

<div class="alert alert-success">
<b>Data Locations:</b><br>
- Raw Data: <code>/srv/data/meg-m2-rnaseq/Data/fastq/raw/</code>
- Genome Annotation: <code>/srv/data/meg-m2-rnaseq/Genomes/Mmu/GRCm39/extracted/mouse genome_annotation-M35.gtf</code>
- Fastp Results: <code>/srv/data/meg-m2-rnaseq/Results/fastp/</code>
</div>"""

# Environment setup
env_md = """## 1. Environment Setup

<div class="alert alert-info">
The required tools are pre-installed in the Plasmabio environment:
- FastQC (v0.12.1) - Quality control
- MultiQC (v1.13) - Aggregate reports
- fastp (v0.23.1) - Read preprocessing
- STAR (v2.7.11a) - Read mapping
- samtools (v1.18) - BAM file handling
</div>"""

# Create directories
setup_code = """# Cell 1: Create working directories
mkdir -p ~/rnaseq/results/{fastqc,star}
cd ~/rnaseq"""

# Quality assessment
qc_md = """## 2. Raw Data Quality Assessment

<div class="alert alert-info">
We'll first examine the quality of raw sequencing data using FastQC.
We'll start by analyzing the first two samples, then provide commands for processing all samples.
</div>"""

fastqc_code = """# Cell 2: Run FastQC on first two samples
cd ~/rnaseq
fastqc -o results/fastqc -t 10 \\
  $(ls /srv/data/meg-m2-rnaseq/Data/fastq/raw/*_R{1,2}.fastq.gz | head -n 4)"""

fastqc_all_code = """# (Cell 3): Command to run FastQC on all samples
cd ~/rnaseq
fastqc -o results/fastqc -t 10 \\
  /srv/data/meg-m2-rnaseq/Data/fastq/raw/*.fastq.gz"""

# Read preprocessing
preproc_md = """## 3. Read Preprocessing

<div class="alert alert-info">
We use fastp to:
- Trim low quality bases
- Remove adapter sequences
- Filter out poor quality reads

<b>Note:</b> Pre-processed results for all samples are available in:
<code>/srv/data/meg-m2-rnaseq/Results/fastp/</code>
</div>"""

fastp_code = """# Cell 4: Process first two samples with fastp
cd ~/rnaseq
for sample in $(ls /srv/data/meg-m2-rnaseq/Data/fastq/raw/*_R1.fastq.gz | head -n 2); do
    base=$(basename $sample _R1.fastq.gz)
    echo "Processing $base..."
    fastp \\
        -i ${sample} \\
        -I ${sample%_R1.fastq.gz}_R2.fastq.gz \\
        -o /srv/data/meg-m2-rnaseq/Results/fastp/${base}_R1_cleaned.fastq.gz \\
        -O /srv/data/meg-m2-rnaseq/Results/fastp/${base}_R2_cleaned.fastq.gz \\
        --html /srv/data/meg-m2-rnaseq/Results/fastp/${base}_report.html \\
        --json /srv/data/meg-m2-rnaseq/Results/fastp/${base}_report.json \\
        --thread 10
done"""

fastp_all_code = """# (Cell 5): Command to process all samples with fastp
cd ~/rnaseq
for sample in /srv/data/meg-m2-rnaseq/Data/fastq/raw/*_R1.fastq.gz; do
    base=$(basename $sample _R1.fastq.gz)
    echo "Processing $base..."
    fastp \\
        -i ${sample} \\
        -I ${sample%_R1.fastq.gz}_R2.fastq.gz \\
        -o /srv/data/meg-m2-rnaseq/Results/fastp/${base}_R1_cleaned.fastq.gz \\
        -O /srv/data/meg-m2-rnaseq/Results/fastp/${base}_R2_cleaned.fastq.gz \\
        --html /srv/data/meg-m2-rnaseq/Results/fastp/${base}_report.html \\
        --json /srv/data/meg-m2-rnaseq/Results/fastp/${base}_report.json \\
        --thread 10
done"""

# Post-processing QC
postqc_md = """## 4. Post-processing Quality Check

<div class="alert alert-info">
We'll run FastQC on the cleaned reads and generate a MultiQC report combining all quality metrics.
We'll use the pre-processed fastp results available in: <code>/srv/data/meg-m2-rnaseq/Results/fastp/</code>
</div>"""

postqc_code = """# Cell 6: Run FastQC on cleaned reads (first two samples)
cd ~/rnaseq
fastqc -o results/fastqc -t 10 \\
  $(ls /srv/data/meg-m2-rnaseq/Results/fastp/*_cleaned.fastq.gz | head -n 4)"""

multiqc_code = """# Cell 7: Generate MultiQC report for all samples
cd ~/rnaseq
# Create symbolic link to fastp results
ln -sf /srv/data/meg-m2-rnaseq/Results/fastp .
# Generate report including all samples
multiqc -o results/multiqc \\
  results/fastqc \\
  /srv/data/meg-m2-rnaseq/Results/fastp"""

# Read mapping
mapping_md = """## 5. Read Mapping

<div class="alert alert-info">
We'll map the cleaned reads to the mouse reference genome using STAR.
We'll demonstrate the mapping process on the first two samples.
The genome annotation file is located at:
<code>/srv/data/meg-m2-rnaseq/Genomes/Mmu/GRCm39/extracted/mouse genome_annotation-M35.gtf</code>
</div>"""

genome_prep_code = """# (Cell 8): Download and prepare reference genome
cd ~/rnaseq
mkdir -p reference
cd reference
wget https://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"""

index_code = """# (Cell 9): Index reference genome
cd ~/rnaseq/reference
STAR --runMode genomeGenerate \\
     --genomeDir star_index \\
     --genomeFastaFiles Mus_musculus.GRCm39.dna.primary_assembly.fa \\
     --sjdbGTFfile /srv/data/meg-m2-rnaseq/Genomes/Mmu/GRCm39/extracted/mouse\\ genome_annotation-M35.gtf \\
     --runThreadN 10"""

mapping_code = """# Cell 10: Map reads for first two samples
cd ~/rnaseq
for sample in $(ls /srv/data/meg-m2-rnaseq/Results/fastp/*_R1_cleaned.fastq.gz | head -n 2); do
    base=$(basename $sample _R1_cleaned.fastq.gz)
    echo "Mapping $base..."
    STAR --genomeDir reference/star_index \\
         --readFilesIn ${sample} ${sample%_R1_cleaned.fastq.gz}_R2_cleaned.fastq.gz \\
         --readFilesCommand zcat \\
         --outFileNamePrefix results/star/${base}_ \\
         --outSAMtype BAM SortedByCoordinate \\
         --runThreadN 10
done"""

mapping_all_code = """# (Cell 11): Command to map all samples
# Note: For the course, we'll provide pre-mapped BAM files for the remaining samples
cd ~/rnaseq
for sample in /srv/data/meg-m2-rnaseq/Results/fastp/*_R1_cleaned.fastq.gz; do
    base=$(basename $sample _R1_cleaned.fastq.gz)
    echo "Mapping $base..."
    STAR --genomeDir reference/star_index \\
         --readFilesIn ${sample} ${sample%_R1_cleaned.fastq.gz}_R2_cleaned.fastq.gz \\
         --readFilesCommand zcat \\
         --outFileNamePrefix results/star/${base}_ \\
         --outSAMtype BAM SortedByCoordinate \\
         --runThreadN 10
done"""

index_bam_code = """# Cell 12: Index BAM files
cd ~/rnaseq
for bam in results/star/*_Aligned.sortedByCoord.out.bam; do
    samtools index -@ 10 $bam
done"""

# Add cells to notebook
cells = [
    nbf.v4.new_markdown_cell(title_md),
    nbf.v4.new_markdown_cell(env_md),
    nbf.v4.new_code_cell(setup_code),
    nbf.v4.new_markdown_cell(qc_md),
    nbf.v4.new_code_cell(fastqc_code),
    nbf.v4.new_raw_cell(fastqc_all_code),
    nbf.v4.new_markdown_cell(preproc_md),
    nbf.v4.new_code_cell(fastp_code),
    nbf.v4.new_raw_cell(fastp_all_code),
    nbf.v4.new_markdown_cell(postqc_md),
    nbf.v4.new_code_cell(postqc_code),
    nbf.v4.new_code_cell(multiqc_code),
    nbf.v4.new_markdown_cell(mapping_md),
    nbf.v4.new_raw_cell(genome_prep_code),
    nbf.v4.new_raw_cell(index_code),
    nbf.v4.new_code_cell(mapping_code),
    nbf.v4.new_raw_cell(mapping_all_code),
    nbf.v4.new_code_cell(index_bam_code)
]

nb.cells = cells

# Set kernel spec
nb.metadata.kernelspec = {
    "display_name": "Bash",
    "language": "bash",
    "name": "bash"
}

# Write the notebook
with open('simplified_notebook.ipynb', 'w') as f:
    nbf.write(nb, f)
