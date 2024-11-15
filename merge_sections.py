import nbformat
import re

def extract_sections(notebook_file, sections):
    """Extract specified sections from notebook"""
    print(f"\nAnalyzing notebook: {notebook_file}")
    with open(notebook_file, 'r', encoding='utf-8') as f:
        nb = nbformat.read(f, as_version=4)

    extracted_cells = []
    current_section = None
    section_content = []
    in_section = False

    print("\nSearching for sections:")
    for section in sections:
        print(f"- {section}")

    print("\nFound markdown headers:")
    for cell in nb.cells:
        if cell.cell_type == 'markdown':
            if re.match(r'^#+\s+', cell.source):
                print(f"Header found: {cell.source.splitlines()[0]}")

    print("\nProcessing sections...")
    for cell in nb.cells:
        if cell.cell_type == 'markdown':
            found_new_section = False
            # Check if this cell starts a new section
            for section in sections:
                # More flexible section matching
                section_pattern = section.replace(" - ", "[ -]+").replace(".", r"\.").replace("featuresCounts", r"(?:<code>)?featuresCounts(?:</code>)?")
                if (section in cell.source or
                    re.search(rf"#+ *{re.escape(section)}", cell.source, re.IGNORECASE) or
                    re.search(rf"#+ *{section_pattern}", cell.source, re.IGNORECASE) or
                    # Handle HTML-formatted headers
                    re.search(rf"#+ *{section.replace('featuresCounts', '<code>featuresCounts</code>')}", cell.source, re.IGNORECASE)):
                    print(f"\nFound section: {section}")
                    if current_section and section_content:
                        extracted_cells.extend(section_content)
                    current_section = section
                    section_content = [cell]
                    in_section = True
                    found_new_section = True
                    break

            # Only check for section end if we didn't just start a new section
            if not found_new_section and in_section and re.match(r'^#+\s+\d+[\.-]', cell.source):
                if current_section and section_content:
                    print(f"\nEnding section: {current_section}")
                    extracted_cells.extend(section_content)
                    current_section = None
                    section_content = []
                    in_section = False
        elif in_section:
            section_content.append(cell)

    # Add the last section if it exists
    if current_section and section_content:
        print(f"\nAdding final section: {current_section}")
        extracted_cells.extend(section_content)

    print(f"\nTotal cells extracted: {len(extracted_cells)}")
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

    if not pipe06_cells:
        print("Warning: No sections were extracted. Check section names and source file.")
        return

    # Add a markdown cell to separate the merged content
    separator_cell = nbformat.v4.new_markdown_cell(source="\n## Merged Sections from Pipe_06\n")
    target_nb.cells.append(separator_cell)

    # Add extracted cells
    target_nb.cells.extend(pipe06_cells)

    # Save modified notebook
    with open(target_file, 'w', encoding='utf-8') as f:
        nbformat.write(target_nb, f)

    print("Sections merged successfully")

if __name__ == "__main__":
    merge_sections("PS5-2024-merged.ipynb", "Pipe_06-bash_reads-counts-pseudomapping.ipynb")
