import nbformat

def create_markdown_cell(content):
    return nbformat.v4.new_markdown_cell(content)

def create_code_cell(content):
    return nbformat.v4.new_code_cell(content)

def add_documentation_section():
    return create_markdown_cell("""## Documentation and References

### Tool Documentation
- [featureCounts Documentation](https://subread.sourceforge.net/featureCounts.html)
- [Salmon Documentation](https://salmon.readthedocs.io/en/latest/)

### Parameter References
- featureCounts parameters: [Manual](http://bioinf.wehi.edu.au/featureCounts/)
- Salmon parameters: [Parameters](https://salmon.readthedocs.io/en/latest/salmon.html#parameters)""")

def add_troubleshooting_section():
    return create_markdown_cell("""## Troubleshooting Guide

Common Issues and Solutions:

1. **Insufficient Memory**
   - Symptom: Process killed or memory error
   - Solution: Reduce number of threads or batch size

2. **Missing Input Files**
   - Symptom: File not found error
   - Solution: Verify file paths and permissions

3. **Invalid File Format**
   - Symptom: Parse error in featureCounts/Salmon
   - Solution: Check input file format and integrity

4. **Resource Issues**
   - Symptom: Slow processing or timeouts
   - Solution: Adjust resource allocation""")

def add_error_handling(code):
    """Add error handling to code cells"""
    return f"""try:
    # Progress indicator
    print("Starting process...")

    {code}

    # Success indicator
    print("Process completed successfully!")
except FileNotFoundError as e:
    print(f"Error: Input file not found - {{e}}")
    raise
except Exception as e:
    print(f"Error occurred: {{e}}")
    raise"""

def add_missing_improvements(notebook_file):
    print(f"Adding missing improvements to: {notebook_file}")

    # Load notebook
    with open(notebook_file, 'r', encoding='utf-8') as f:
        nb = nbformat.read(f, as_version=4)

    # Add documentation section after first markdown cell
    for i, cell in enumerate(nb.cells):
        if cell.cell_type == 'markdown':
            nb.cells.insert(i + 1, add_documentation_section())
            break

    # Add troubleshooting section before the last cell
    nb.cells.insert(len(nb.cells) - 1, add_troubleshooting_section())

    # Add error handling and progress indicators to code cells
    for cell in nb.cells:
        if cell.cell_type == 'code' and not cell.source.startswith('try:'):
            cell.source = add_error_handling(cell.source)

    # Save modified notebook
    with open(notebook_file, 'w', encoding='utf-8') as f:
        nbformat.write(nb, f)

    print("Missing improvements added successfully")

if __name__ == "__main__":
    add_missing_improvements("PS5-2024-merged.ipynb")
