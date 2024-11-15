import nbformat
import sys
from pathlib import Path

def analyze_notebook(notebook_path, check_sections=None):
    """Analyze a Jupyter notebook and print its structure."""
    print(f"\nAnalyzing notebook: {notebook_path}")

    with open(notebook_path, 'r', encoding='utf-8') as f:
        nb = nbformat.read(f, as_version=4)

    print("\nStructure:")
    print("-" * 50)

    cell_count = {'markdown': 0, 'code': 0, 'raw': 0}
    sections = []
    current_section = None
    code_cells_content = []
    found_sections = set()
    current_full_section = []

    for idx, cell in enumerate(nb.cells):
        cell_count[cell.cell_type] += 1

        if cell.cell_type == 'markdown':
            lines = cell.source.split('\n')
            for line in lines:
                if line.startswith('#'):
                    section_text = line.strip().replace('**', '').replace('<code>', '').replace('</code>', '')
                    if check_sections:
                        for required in check_sections:
                            # Remove formatting and do case-insensitive comparison
                            clean_required = required.replace('**', '').lower()
                            clean_section = section_text.lower()
                            if clean_required in clean_section:
                                found_sections.add(required)
                    sections.append(section_text)
        elif cell.cell_type == 'code':
            # Check code cell numbering format
            if not cell.source.strip().startswith('## Code cell'):
                print(f"\nWARNING: Code cell {idx} does not follow numbering format!")
                print(f"Content preview: {cell.source[:50]}...")
            code_cells_content.append(f"Cell {idx}: {cell.source[:100]}...")

    print("\nCell Statistics:")
    print(f"Total cells: {len(nb.cells)}")
    for cell_type, count in cell_count.items():
        print(f"{cell_type.capitalize()} cells: {count}")

    print("\nMain Sections:")
    for section in sections:
        print(f"- {section}")

    if check_sections:
        print("\nRequired Sections Check:")
        for section in check_sections:
            status = "✓ Found" if section in found_sections else "✗ Missing"
            print(f"{status}: {section}")

    print("\nCode Cell Previews:")
    for preview in code_cells_content:
        print(preview)

if __name__ == "__main__":
    required_sections = [
        "1.3.2 - Running featuresCounts on multiple samples",
        "2 - Pseudo-mapping with Salmon"
    ]

    merged_path = "/home/ubuntu/workspace/rnaseq/RNAseq_Plasma/PS5-2024-merged.ipynb"
    print("=" * 80)
    print("Analyzing Merged Notebook")
    print("=" * 80)
    analyze_notebook(merged_path, check_sections=required_sections)
