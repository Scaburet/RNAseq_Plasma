import nbformat
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

def add_progressive_loop(notebook_path):
    # Read the existing notebook
    with open(notebook_path, 'r', encoding='utf-8') as f:
        nb = nbformat.read(f, as_version=4)

    # Add progressive loop building section
    loop_cells = [
        new_markdown_cell("""
## Progressive Loop Building for Qualimap Analysis

In this section, we'll build a loop to run Qualimap on all samples step by step. This approach helps understand how to automate repetitive tasks in bash.

### Step 1: Understanding the Basic Command
First, let's look at the basic Qualimap command we used for a single sample:
"""),

        new_code_cell("""# Example for a single sample
qualimap bamqc \
    -bam Results/samtools/sample1_sorted.bam \
    -outdir Results/qualimap/sample1 \
    --java-mem-size=6G"""),

        new_markdown_cell("""
### Step 2: Identifying Variable Parts
In the command above, we can identify two main parts that change for each sample:
1. Input BAM file path (`Results/samtools/sample1_sorted.bam`)
2. Output directory path (`Results/qualimap/sample1`)

The sample name is the key variable that changes in both paths.
"""),

        new_markdown_cell("""
### Step 3: Creating a List of Sample Names
First, we'll create a list of our sample names. We can do this by listing the sorted BAM files and extracting the sample names:
"""),

        new_code_cell("""# List sample names
ls Results/samtools/*_sorted.bam | sed 's|Results/samtools/||' | sed 's|_sorted.bam||'"""),

        new_markdown_cell("""
### Step 4: Building the Loop - Basic Structure
Now, let's create a simple loop structure that will iterate over our samples:
"""),

        new_code_cell("""# Basic loop structure
for sample in $(ls Results/samtools/*_sorted.bam | sed 's|Results/samtools/||' | sed 's|_sorted.bam||')
do
    echo "Processing sample: $sample"
done"""),

        new_markdown_cell("""
### Step 5: Adding the Qualimap Command
Now we'll add the Qualimap command inside our loop, using variables for the sample-specific parts:
"""),

        new_code_cell("""# Complete loop with Qualimap command
for sample in $(ls Results/samtools/*_sorted.bam | sed 's|Results/samtools/||' | sed 's|_sorted.bam||')
do
    echo "Processing sample: $sample"

    # Create output directory if it doesn't exist
    mkdir -p Results/qualimap/${sample}

    # Run Qualimap
    qualimap bamqc \
        -bam Results/samtools/${sample}_sorted.bam \
        -outdir Results/qualimap/${sample} \
        --java-mem-size=6G
done"""),

        new_markdown_cell("""
### Important Notes:
1. The `mkdir -p` command ensures our output directory exists
2. We use `${sample}` to clearly separate the variable name in the paths
3. The `--java-mem-size=6G` parameter is important due to the memory constraints (6GB per student)
4. Using `echo` statements helps track progress during execution

### Exercise for Students:
Try modifying the loop to:
1. Add error checking
2. Include progress percentage
3. Save a log file of the analysis
""")
    ]

    # Add the new cells to the notebook
    nb.cells.extend(loop_cells)

    # Write the updated notebook
    with open(notebook_path, 'w', encoding='utf-8') as f:
        nbformat.write(nb, f)

    print("Successfully added progressive loop building section to the notebook")

if __name__ == "__main__":
    notebook_path = "/home/ubuntu/workspace/rnaseq/CEA_rnaseq/temp_repo/PS4-2024-merged.ipynb"
    add_progressive_loop(notebook_path)
