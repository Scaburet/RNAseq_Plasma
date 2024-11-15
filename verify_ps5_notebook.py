import nbformat
import re

def verify_notebook(filename):
    print(f"Verifying notebook: {filename}")

    # Check filename
    if not filename.startswith("PS5-2024"):
        print("❌ Filename does not start with PS5-2024")
        return False
    print("✅ Filename starts with PS5-2024")

    # Load notebook
    with open(filename, 'r', encoding='utf-8') as f:
        nb = nbformat.read(f, as_version=4)

    # Initialize verification flags
    found_featurecounts_section = False
    found_salmon_section = False
    code_cell_count = 0
    all_cells_numbered = True

    # Check cell contents and numbering
    for cell in nb.cells:
        # Check for required sections
        if cell.cell_type == "markdown":
            if "1.3.2" in cell.source and "Running featureCounts on Multiple Samples" in cell.source:
                found_featurecounts_section = True
            elif "2" in cell.source and "Pseudo-mapping with Salmon" in cell.source:
                found_salmon_section = True

        # Check code cell numbering
        if cell.cell_type == "code":
            code_cell_count += 1
            if not re.match(r'^## Code cell \d+ ##', cell.source):
                print(f"❌ Code cell {code_cell_count} is not properly numbered")
                all_cells_numbered = False

    # Print verification results
    print("\nVerification Results:")
    print(f"✅ Found featureCounts section: {found_featurecounts_section}")
    print(f"✅ Found Salmon section: {found_salmon_section}")
    print(f"✅ Total code cells: {code_cell_count}")
    print(f"{'✅' if all_cells_numbered else '❌'} All code cells properly numbered")

    # Overall verification
    success = (found_featurecounts_section and
              found_salmon_section and
              all_cells_numbered)

    print(f"\nOverall verification: {'✅ PASSED' if success else '❌ FAILED'}")
    return success

if __name__ == "__main__":
    verify_notebook("PS5-2024-merged.ipynb")
