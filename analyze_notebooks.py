import nbformat
import json
from pathlib import Path

def analyze_notebook(notebook_path):
    """Analyze a Jupyter notebook and return its structure and content."""
    try:
        with open(notebook_path, 'r', encoding='utf-8') as f:
            nb = nbformat.read(f, as_version=4)

        cell_types = {'markdown': 0, 'code': 0, 'raw': 0}
        content_summary = []

        for cell in nb.cells:
            cell_types[cell.cell_type] += 1
            if cell.cell_type == 'markdown':
                # Get first line of markdown cell to understand section
                first_line = cell.source.split('\n')[0][:100]
                content_summary.append(f"MD: {first_line}")
            elif cell.cell_type == 'code':
                # Get first line of code cell
                first_line = cell.source.split('\n')[0][:100]
                content_summary.append(f"Code: {first_line}")
            elif cell.cell_type == 'raw':
                first_line = cell.source.split('\n')[0][:100]
                content_summary.append(f"Raw: {first_line}")

        print(f"\nAnalysis for {Path(notebook_path).name}:")
        print("Cell count by type:", json.dumps(cell_types, indent=2))
        print("\nContent Structure:")
        for item in content_summary:
            print(f"- {item}")

    except Exception as e:
        print(f"Error analyzing {notebook_path}: {str(e)}")

# Analyze PS4 notebook
ps4_path = "/home/ubuntu/attachments/PS4-mappingOutput-bash-2022-run.ipynb"
analyze_notebook(ps4_path)

# Analyze Pipe_05 notebook (using the found path)
pipe05_path = "/home/ubuntu/workspace/rnaseq/Pipe_05-bash_mapping-quality.ipynb"
analyze_notebook(pipe05_path)
