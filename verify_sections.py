import nbformat
import re

def verify_sections(notebook_file):
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
                # Check for exact match or similar section header
                if section in cell.source or re.search(rf"#+ *{re.escape(section)}", cell.source, re.IGNORECASE):
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

if __name__ == "__main__":
    verify_sections("PS5-2024-merged.ipynb")
