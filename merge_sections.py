import nbformat
import re

def extract_sections(notebook_file, sections):
    """Extract specified sections from notebook"""
    with open(notebook_file, 'r', encoding='utf-8') as f:
        nb = nbformat.read(f, as_version=4)

    extracted_cells = []
    current_section = None
    section_content = []

    for cell in nb.cells:
        if cell.cell_type == 'markdown':
            # Check if this cell starts a new section
            for section in sections:
                if section in cell.source or re.search(rf"#+ *{re.escape(section)}", cell.source, re.IGNORECASE):
                    if current_section and section_content:
                        extracted_cells.extend(section_content)
                    current_section = section
                    section_content = [cell]
                    break
        elif current_section:
            section_content.append(cell)

    # Add the last section if it exists
    if section_content:
        extracted_cells.extend(section_content)

    return extracted_cells

def merge_sections(target_file, pipe06_file):
    """Merge required sections into target notebook"""
    print(f"Merging sections into: {target_file}")

    required_sections = [
        "1.3.2 - Running featuresCounts on multiple samples",
        "2 - Pseudo-mapping with Salmon"
    ]

    # Load target notebook
    with open(target_file, 'r', encoding='utf-8') as f:
        target_nb = nbformat.read(f, as_version=4)

    # Extract sections from Pipe_06
    pipe06_cells = extract_sections(pipe06_file, required_sections)

    # Find insertion point (end of notebook)
    target_nb.cells.extend(pipe06_cells)

    # Save modified notebook
    with open(target_file, 'w', encoding='utf-8') as f:
        nbformat.write(target_nb, f)

    print("Sections merged successfully")

if __name__ == "__main__":
    merge_sections("PS5-2024-merged.ipynb", "Pipe_06-bash_reads-counts-pseudomapping.ipynb")
