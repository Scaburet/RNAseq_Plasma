import nbformat
import json

def analyze_notebook_content(notebook_path):
    """Analyze the content and structure of a Jupyter notebook."""
    with open(notebook_path, 'r', encoding='utf-8') as f:
        nb = nbformat.read(f, as_version=4)

    # Count cells by type
    cell_counts = {'markdown': 0, 'code': 0, 'raw': 0}
    sections = []
    current_section = None

    for cell in nb.cells:
        cell_counts[cell.cell_type] += 1

        if cell.cell_type == 'markdown':
            # Check for section headers
            lines = cell.source.split('\n')
            for line in lines:
                if line.startswith('#'):
                    sections.append(line.strip())
                if 'alert' in line:
                    if current_section:
                        current_section['alerts'] = current_section.get('alerts', 0) + 1

    print(f"\nAnalysis for {notebook_path}:")
    print("Cell count by type:", json.dumps(cell_counts, indent=2))
    print("\nMain sections:")
    for section in sections:
        print(f"- {section}")

    # Verify key components
    required_components = [
        "environment setup",
        "qualimap",
        "multiqc",
        "bam",
        "results",
        "genome"
    ]

    content = ' '.join(cell.source.lower() for cell in nb.cells)
    missing_components = [comp for comp in required_components if comp not in content]

    if missing_components:
        print("\nMissing key components:", missing_components)
    else:
        print("\nAll key components present")

# Analyze the merged notebook
analyze_notebook_content('PS4-2024-merged.ipynb')
