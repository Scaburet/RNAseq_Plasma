import nbformat
import re

def verify_sections(notebook_file):
    """Verify required sections are present in notebook"""
    print(f"Verifying required sections in: {notebook_file}")

    required_sections = [
        "1.3.2 - Running featuresCounts on multiple samples",
        "2 - Pseudo-mapping with Salmon"
    ]

    found_sections = {section: False for section in required_sections}

    # Load notebook
    with open(notebook_file, 'r', encoding='utf-8') as f:
        nb = nbformat.read(f, as_version=4)

    # Check each cell for section headers
    for cell in nb.cells:
        if cell.cell_type == 'markdown':
            for section in required_sections:
                # Create HTML-aware pattern
                section_pattern = section.replace(" - ", "[ -]+").replace(".", r"\.").replace("featuresCounts", r"(?:<code>)?featuresCounts(?:</code>)?")
                if (section in cell.source or
                    re.search(rf"#+ *{re.escape(section)}", cell.source, re.IGNORECASE) or
                    re.search(rf"#+ *{section_pattern}", cell.source, re.IGNORECASE) or
                    re.search(rf"#+ *{section.replace('featuresCounts', '<code>featuresCounts</code>')}", cell.source, re.IGNORECASE)):
                    found_sections[section] = True

    # Print results
    print("\nSection Verification Results:")
    all_found = True
    for section, found in found_sections.items():
        status = "✅" if found else "❌"
        print(f"{status} {section}")
        if not found:
            all_found = False

    print(f"\nOverall section verification: {'✅ PASSED' if all_found else '❌ FAILED'}")
    return all_found

def number_code_cells(notebook_file):
    """Add cell numbers to code cells in notebook"""
    print(f"\nNumbering code cells in: {notebook_file}")

    with open(notebook_file, 'r', encoding='utf-8') as f:
        nb = nbformat.read(f, as_version=4)

    code_cell_count = 1
    for cell in nb.cells:
        if cell.cell_type == 'code':
            # Add or update cell number comment at the start
            if not cell.source.startswith('## Code cell'):
                cell.source = f"## Code cell {code_cell_count} ##\n{cell.source}"
            else:
                cell.source = re.sub(r'^## Code cell \d+ ##', f'## Code cell {code_cell_count} ##', cell.source)
            code_cell_count += 1

    with open(notebook_file, 'w', encoding='utf-8') as f:
        nbformat.write(nb, f)

    print(f"Numbered {code_cell_count-1} code cells")
    return code_cell_count-1

if __name__ == "__main__":
    notebook_file = "PS5-2024-merged.ipynb"
    verify_sections(notebook_file)
    number_code_cells(notebook_file)
