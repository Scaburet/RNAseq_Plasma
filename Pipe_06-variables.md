# Minimum Required Variables for Pipe_06

## Core System Variables
```bash
# Base directories
WORK_DIR="/srv/home/scaburet/meg_m2_rnaseq_bash"
DATA_DIR="/srv/data/meg-m2-rnaseq"

# Input/Output directories
RESULTS_DIR="${WORK_DIR}/Results"
BAM_DIR="${RESULTS_DIR}/star"
COUNTS_DIR="${RESULTS_DIR}/counts"
SALMON_DIR="${RESULTS_DIR}/salmon"

# Reference files
GTF_FILE="${DATA_DIR}/genome/Homo_sapiens.GRCh38.106.gtf"
TRANSCRIPTOME="${DATA_DIR}/genome/Homo_sapiens.GRCh38.cdna.all.fa"
```

## Resource Management
```bash
# System resources
CPU_CORES=4
MAX_RAM="20G"
```

## Tool Parameters
```bash
# featureCounts settings
STRAND_SPECIFIC=0  # 0=unstranded, 1=stranded, 2=reversely stranded
MIN_MAPPING_QUALITY=10

# Salmon settings
LIBTYPE="A"  # Automatic detection of library type
```

## Notes
1. All paths are absolute and Plasmabio-specific
2. Resource limits set for shared environment
3. Tool parameters optimized for educational use
4. Directory structure simplified for clarity
