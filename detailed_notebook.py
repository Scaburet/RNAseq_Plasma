#!/usr/bin/env python
# coding: utf-8

# # RNA-seq Data Analysis Pipeline
# ## Date: Tuesday 26/11/2024
# 
# This notebook contains the complete workflow for RNA-seq data analysis, including quality control, preprocessing, and mapping steps. The notebook is designed to be self-explanatory and includes detailed documentation for future reference.
# 
# ## Overview
# This pipeline processes paired-end RNA-seq data from mouse samples. We will:
# 1. Check the quality of raw reads using FastQC
# 2. Preprocess reads using fastp
# 3. Verify the quality of processed reads
# 4. Download and index the reference genome
# 5. Map reads to the reference genome using STAR
# 6. Generate comprehensive quality reports with MultiQC
# 
# ## Data Organization
# - Raw data location: `/srv/data/meg-m2-rnaseq/Data/fastq/raw/`
# - Genome annotation: `/srv/data/meg-m2-rnaseq/Genomes/Mmu/GRCm39/extracted/`
# - Results will be organized in tool-specific subfolders:
#   - FastQC results: `Results/fastqc/`
#   - fastp results: `Results/fastp/`
#   - STAR mapping results: `Results/star/`
#   - MultiQC results: `Results/multiqc/`
# 
# ## Required Software
# - FastQC v0.12.1
# - MultiQC v1.13
# - fastp v0.23.1
# - STAR v2.7.11a
# - samtools v1.18
# 
# ## Important Notes
# - This notebook demonstrates the workflow on the first two samples
# - Commands for processing all samples are provided as commented code
# - Results are organized in tool-specific directories
# - The pipeline is designed for mouse RNA-seq data
# 

# In[1]:


get_ipython().run_cell_magic('bash', '', '# Create Results directory structure\nmkdir -p Results/{fastqc,fastp,star,multiqc}\n\n# Display the created directory structure\ntree Results/\n')


# ## Reference Genome Download and Indexing
# 
# ### Overview
# The mouse reference genome and its annotation are fundamental for RNA-seq analysis:
# - Reference genome: DNA sequence assembly (GRCm39)
# - Gene annotation: Gene models and features (GTF format)
# - Index: Optimized data structure for rapid read mapping
# 
# ### Understanding Reference Files
# 1. **Reference Genome (FASTA)**
#    - Contains chromosome sequences
#    - Represents mouse genome assembly GRCm39
#    - Downloaded from Ensembl
#    - Essential for read alignment
# 
# 2. **Gene Annotation (GTF)**
#    - Defines gene structures
#    - Contains exon/intron boundaries
#    - Identifies gene features
#    - Required for splice-aware mapping
# 
# 3. **STAR Index**
#    - Optimized search structure
#    - Enables fast read mapping
#    - Contains splice junction information
#    - Memory-intensive to generate
# 
# ### Resource Requirements
# - Disk Space: ~8GB for genome and index
# - Memory: ~32GB recommended for indexing
# - CPU: Multi-threaded process (using 10 threads)
# - Time: 30-60 minutes depending on resources
# 
# ### Key Parameters for STAR Indexing
# - `--runMode genomeGenerate`: Index generation mode
# - `--genomeDir`: Index output directory
# - `--genomeFastaFiles`: Reference genome path
# - `--sjdbGTFfile`: Gene annotation file
# - `--runThreadN`: Number of threads
# - `--sjdbOverhang`: Read length - 1 (default: 100)
# 
# ### Documentation Resources
# - [Ensembl Mouse Genome](https://www.ensembl.org/Mus_musculus/Info/Index)
# - [STAR Manual - Genome Generation](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf#page=7)
# - [GTF Format Specification](https://www.ensembl.org/info/website/upload/gff.html)
# - [GRCm39 Assembly Info](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.27)
# 

# In[2]:


get_ipython().run_cell_magic('bash', '', '# Download and prepare reference genome\necho "Creating reference directories..."\nmkdir -p reference/star_index\n\necho "Downloading mouse reference genome..."\ncd reference\nwget https://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz\ngunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz\n\necho "Creating STAR index..."\nSTAR --runMode genomeGenerate \\\n     --genomeDir star_index \\\n     --genomeFastaFiles Mus_musculus.GRCm39.dna.primary_assembly.fa \\\n     --sjdbGTFfile "/srv/data/meg-m2-rnaseq/Genomes/Mmu/GRCm39/extracted/mouse genome_annotation-M35.gtf" \\\n     --runThreadN 10 \\\n     --sjdbOverhang 100\n\necho "Checking index files..."\nls -lh star_index/\necho "Index generation complete."\ncd ..\n')


# ## Quality Control with FastQC
# 
# ### Overview
# FastQC performs quality control checks on raw sequence data. Understanding these quality metrics is crucial for:
# - Identifying sequencing problems
# - Detecting contamination
# - Ensuring data reliability
# - Planning preprocessing steps
# 
# ### Quality Metrics Explained
# 1. **Basic Statistics**
#    - Total sequences
#    - Sequence length
#    - GC content percentage
#    - Important for sample validation
# 
# 2. **Per Base Sequence Quality**
#    - Quality scores across all bases
#    - Look for declining quality at read ends
#    - Phred scores > 28 are considered good
#    - Common drop in quality at read ends
# 
# 3. **Per Sequence Quality Scores**
#    - Overall quality distribution
#    - Most sequences should have mean quality > 27
#    - Identifies problematic sequences
# 
# 4. **Sequence Duplication Levels**
#    - PCR duplication assessment
#    - Expected higher in RNA-seq
#    - Helps identify library complexity
# 
# 5. **Overrepresented Sequences**
#    - Adapter contamination check
#    - Highly expressed genes
#    - rRNA contamination
# 
# ### Important Parameters
# - `-o`: Output directory for reports
# - `-t`: Number of threads (we use 10)
# - `--noextract`: Keep reports zipped
# - `--quiet`: Suppress progress messages
# 
# ### Documentation Resources
# - [FastQC Documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
# - [Example Reports](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html)
# - [Quality Score Encoding](https://en.wikipedia.org/wiki/FASTQ_format#Quality)
# - [RNA-seq QC Best Practices](https://rna-seqblog.com/quality-control-tools/)
# 

# In[3]:


get_ipython().run_cell_magic('bash', '', '# Run FastQC on first two samples\necho "Creating FastQC output directory..."\nmkdir -p Results/fastqc\n\nDATA_DIR="/srv/data/meg-m2-rnaseq/Data/fastq/raw"\necho "Running FastQC on first two samples..."\nfastqc -o Results/fastqc -t 10 \\\n    ${DATA_DIR}/sample1_R1.fastq.gz \\\n    ${DATA_DIR}/sample1_R2.fastq.gz \\\n    ${DATA_DIR}/sample2_R1.fastq.gz \\\n    ${DATA_DIR}/sample2_R2.fastq.gz\n\necho "FastQC analysis complete. Checking results..."\nls -lh Results/fastqc/\n\n# Command for all samples (shown but not executed)\n: \'\necho "Running FastQC on all samples..."\nfastqc -o Results/fastqc -t 10 ${DATA_DIR}/*.fastq.gz\n\'\n')


# ## Read Preprocessing with fastp
# 
# ### Overview
# fastp performs quality control and preprocessing of sequencing reads. This step is crucial for:
# - Removing technical artifacts
# - Improving mapping accuracy
# - Reducing false positives in downstream analysis
# - Generating detailed QC reports
# 
# ### Understanding the Preprocessing Steps
# 1. **Adapter Trimming**
#    - Removes sequencing adapters automatically
#    - Prevents misalignments
#    - Improves mapping efficiency
#    - Uses overlap analysis for adapter detection
# 
# 2. **Quality Filtering**
#    - Removes low-quality bases from read ends
#    - Filters out poor quality reads
#    - Uses sliding window approach
#    - Configurable quality thresholds
# 
# 3. **Base Correction (Paired-End Only)**
#    - Uses overlapping information
#    - Corrects mismatched bases
#    - Improves sequence accuracy
#    - Only applies to overlapping regions
# 
# 4. **Length Filtering**
#    - Removes too short reads
#    - Prevents spurious mappings
#    - Configurable minimum length
#    - Important for mapping accuracy
# 
# ### Quality Report Components
# 1. **Before Processing**
#    - Initial quality metrics
#    - Base composition
#    - Quality score distribution
#    - Adapter content
# 
# 2. **After Processing**
#    - Final read statistics
#    - Quality improvements
#    - Filtering statistics
#    - Base quality distribution
# 
# ### Key Parameters Explained
# - `--in1/--in2`: Input paired-end files
# - `--out1/--out2`: Output cleaned files
# - `--html`: HTML report location
# - `--json`: JSON format metrics
# - `--thread`: CPU threads (10 recommended)
# - `--qualified_quality_phred`: Minimum base quality (default: 15)
# - `--length_required`: Minimum read length (default: 36)
# - `--cut_front`: Trim low quality bases from read start
# - `--cut_tail`: Trim low quality bases from read end
# - `--cut_mean_quality`: Average quality threshold
# 
# ### Documentation Resources
# - [fastp GitHub](https://github.com/OpenGene/fastp)
# - [Parameters Guide](https://github.com/OpenGene/fastp#all-options)
# - [Algorithm Details](https://github.com/OpenGene/fastp#algorithm)
# - [Publication](https://doi.org/10.1093/bioinformatics/bty560)
# 

# In[4]:


get_ipython().run_cell_magic('bash', '', '# Process first two samples with fastp\necho "Creating fastp output directory..."\nmkdir -p Results/fastp\n\nDATA_DIR="/srv/data/meg-m2-rnaseq/Data/fastq/raw"\nfor sample in sample1 sample2; do\n    echo "Processing ${sample} with fastp..."\n\n    # Run fastp with detailed parameters\n    fastp --in1 ${DATA_DIR}/${sample}_R1.fastq.gz \\\n          --in2 ${DATA_DIR}/${sample}_R2.fastq.gz \\\n          --out1 Results/fastp/${sample}_R1_cleaned.fastq.gz \\\n          --out2 Results/fastp/${sample}_R2_cleaned.fastq.gz \\\n          --html Results/fastp/${sample}_report.html \\\n          --json Results/fastp/${sample}_report.json \\\n          --thread 10 \\\n          --qualified_quality_phred 15 \\\n          --length_required 36 \\\n          --cut_front \\\n          --cut_tail \\\n          --cut_mean_quality 20\n\n    echo "Checking output files for ${sample}..."\n    ls -lh Results/fastp/${sample}*\n\n    echo "Summary of processing for ${sample}:"\n    echo "--------------------------------"\n    echo "Original files:"\n    ls -lh ${DATA_DIR}/${sample}_R*.fastq.gz\n    echo "Processed files:"\n    ls -lh Results/fastp/${sample}_R*_cleaned.fastq.gz\n    echo "Reports:"\n    ls -lh Results/fastp/${sample}_report.*\n    echo "--------------------------------"\ndone\n\n# Command for all samples (shown but not executed)\n: \'\necho "Processing all samples with fastp..."\nfor sample in sample*; do\n    fastp --in1 ${DATA_DIR}/${sample}_R1.fastq.gz \\\n          --in2 ${DATA_DIR}/${sample}_R2.fastq.gz \\\n          --out1 Results/fastp/${sample}_R1_cleaned.fastq.gz \\\n          --out2 Results/fastp/${sample}_R2_cleaned.fastq.gz \\\n          --html Results/fastp/${sample}_report.html \\\n          --json Results/fastp/${sample}_report.json \\\n          --thread 10 \\\n          --qualified_quality_phred 15 \\\n          --length_required 36 \\\n          --cut_front \\\n          --cut_tail \\\n          --cut_mean_quality 20\ndone\n\'\n')


# ## Read Mapping with STAR
# 
# ### Overview
# STAR (Spliced Transcripts Alignment to a Reference) is designed specifically for RNA-seq data:
# - Handles splice junctions efficiently
# - Supports paired-end reads
# - High accuracy and speed
# - Generates comprehensive mapping statistics
# 
# ### Understanding the Mapping Process
# 1. **Seed Search**
#    - Finds Maximum Mappable Prefixes (MMPs)
#    - Uses uncompressed suffix arrays
#    - Extremely fast and memory-efficient
#    - Identifies potential mapping locations
# 
# 2. **Clustering/Stitching**
#    - Clusters and extends seed matches
#    - Accounts for sequencing errors
#    - Handles splice junctions
#    - Identifies chimeric alignments
# 
# 3. **Scoring and Selection**
#    - Evaluates alignment quality
#    - Considers mismatches and gaps
#    - Handles multi-mappers
#    - Assigns mapping quality scores
# 
# 4. **Output Generation**
#    - Creates sorted BAM files
#    - Includes mapping statistics
#    - Records splice junctions
#    - Generates alignment logs
# 
# ### Key Parameters Explained
# - `--runThreadN`: Number of threads (10 recommended)
# - `--genomeDir`: Path to STAR index
# - `--readFilesIn`: Input FASTQ files (R1 R2 for paired-end)
# - `--readFilesCommand`: Command for reading compressed files
# - `--outFileNamePrefix`: Prefix for output files
# - `--outSAMtype`: Output format (BAM SortedByCoordinate)
# - `--outBAMsortingThreadN`: Threads for BAM sorting
# - `--limitBAMsortRAM`: Memory limit for sorting
# 
# ### Understanding Output Files
# 1. **BAM Files**
#    - Contains aligned reads
#    - Sorted by coordinate
#    - Indexed for quick access
#    - Used for downstream analysis
# 
# 2. **Log Files**
#    - Final.out: Summary statistics
#    - Log.out: Detailed progress
#    - SJ.out.tab: Splice junctions
#    - Log.progress.out: Real-time progress
# 
# ### Documentation Resources
# - [STAR Manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)
# - [RNA-STAR GitHub](https://github.com/alexdobin/STAR)
# - [Mapping Statistics Guide](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/)
# - [BAM Format Specification](https://samtools.github.io/hts-specs/SAMv1.pdf)
# 

# In[5]:


get_ipython().run_cell_magic('bash', '', '# Map first two samples with STAR\necho "Creating STAR output directory..."\nmkdir -p Results/star\n\nfor sample in sample1 sample2; do\n    echo "Mapping ${sample} with STAR..."\n\n    # Create sample-specific output directory\n    mkdir -p Results/star/${sample}\n\n    STAR --runThreadN 10 \\\n         --genomeDir reference/star_index \\\n         --readFilesIn Results/fastp/${sample}_R1_cleaned.fastq.gz \\\n                      Results/fastp/${sample}_R2_cleaned.fastq.gz \\\n         --readFilesCommand zcat \\\n         --outFileNamePrefix Results/star/${sample}/ \\\n         --outSAMtype BAM SortedByCoordinate \\\n         --outBAMsortingThreadN 10 \\\n         --limitBAMsortRAM 31532137230 \\\n         --outReadsUnmapped Fastx\n\n    echo "Indexing BAM file for ${sample}..."\n    samtools index Results/star/${sample}/Aligned.sortedByCoord.out.bam\n\n    echo "Checking mapping results for ${sample}..."\n    ls -lh Results/star/${sample}/\n    echo "Mapping statistics:"\n    cat Results/star/${sample}/Log.final.out\n    echo "--------------------------------"\ndone\n\n# Commands for all samples (shown but not executed)\n: \'\necho "Mapping all samples with STAR..."\nfor sample in sample*; do\n    mkdir -p Results/star/${sample}\n\n    STAR --runThreadN 10 \\\n         --genomeDir reference/star_index \\\n         --readFilesIn Results/fastp/${sample}_R1_cleaned.fastq.gz \\\n                      Results/fastp/${sample}_R2_cleaned.fastq.gz \\\n         --readFilesCommand zcat \\\n         --outFileNamePrefix Results/star/${sample}/ \\\n         --outSAMtype BAM SortedByCoordinate \\\n         --outBAMsortingThreadN 10 \\\n         --limitBAMsortRAM 31532137230 \\\n         --outReadsUnmapped Fastx\n\n    samtools index Results/star/${sample}/Aligned.sortedByCoord.out.bam\ndone\n\'\n')


# ## Quality Report Aggregation with MultiQC
# 
# ### Overview
# MultiQC aggregates quality control reports from multiple bioinformatics tools into a single comprehensive report:
# - Combines results from FastQC, fastp, and STAR
# - Creates interactive visualizations
# - Enables easy sample comparison
# - Identifies systematic issues across samples
# 
# ### Understanding MultiQC Reports
# 1. **General Statistics**
#    - Sample-level metrics
#    - Read counts and mapping rates
#    - Quality score distributions
#    - Duplication rates
# 
# 2. **FastQC Integration**
#    - Sequence quality histograms
#    - GC content distribution
#    - Adapter content
#    - Overrepresented sequences
# 
# 3. **fastp Results**
#    - Filtering statistics
#    - Quality improvements
#    - Base quality distributions
#    - Read length distribution
# 
# 4. **STAR Alignment Metrics**
#    - Mapping rates
#    - Uniquely mapped reads
#    - Multi-mapped reads
#    - Splice junction statistics
# 
# ### Report Navigation
# 1. **Interactive Features**
#    - Sortable tables
#    - Zoomable plots
#    - Exportable data
#    - Customizable views
# 
# 2. **Data Interpretation**
#    - Sample comparison
#    - Quality threshold indicators
#    - Warning flags
#    - Failure notifications
# 
# ### Key Parameters
# - `-o`: Output directory location
# - `-f`: Force overwrite existing reports
# - `-d`: Print directory creation
# - `--interactive`: Enable interactive plots
# - `--zip-data-dir`: Compress data directory
# 
# ### Documentation Resources
# - [MultiQC Documentation](https://multiqc.info/)
# - [Available Modules](https://multiqc.info/docs/#available-modules)
# - [Report Examples](https://multiqc.info/examples/rna-seq/multiqc_report.html)
# - [Configuration Guide](https://multiqc.info/docs/#configuring-multiqc)
# 

# In[6]:


get_ipython().run_cell_magic('bash', '', '# Run MultiQC on all results\necho "Creating MultiQC output directory..."\nmkdir -p Results/multiqc\n\necho "Running MultiQC..."\nmultiqc --force \\\n        --interactive \\\n        --verbose \\\n        --dirs \\\n        --dir-pattern \'*\' \\\n        -o Results/multiqc \\\n        Results/fastqc Results/fastp Results/star\n\necho "MultiQC report generation complete."\necho "Report location:"\nls -lh Results/multiqc/multiqc_report.html\n\necho "Available quality reports:"\necho "------------------------"\necho "FastQC reports:"\nls -l Results/fastqc/*.html\necho "------------------------"\necho "fastp reports:"\nls -l Results/fastp/*.html\necho "------------------------"\necho "STAR logs:"\nls -l Results/star/*/Log.final.out\necho "------------------------"\n')

