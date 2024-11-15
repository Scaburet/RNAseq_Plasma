import nbformat
import re

def verify_plasmabio_compatibility(filename):
    print(f"Verifying Plasmabio compatibility for: {filename}")

    # Load notebook
    with open(filename, 'r', encoding='utf-8') as f:
        nb = nbformat.read(f, as_version=4)

    # Patterns to check for IFB references
    ifb_patterns = [
        r'.*\bIFB\b.*',
        r'.*\bifb\b.*',
        r'.*\bifb-core\b.*',
        r'.*\bifb\..*',
        r'/shared/.*',
        r'/ifb/.*',
        r'.*\.ifb\..*'
    ]

    # Patterns to verify Plasmabio paths
    plasmabio_paths = [
        r'/srv/home/\${USER}/meg_m2_rnaseq_bash',
        r'/srv/data/meg-m2-rnaseq'
    ]

    # Compile patterns
    ifb_patterns = [re.compile(pattern, re.IGNORECASE) for pattern in ifb_patterns]
    plasmabio_patterns = [re.compile(pattern) for pattern in plasmabio_paths]

    # Initialize verification results
    ifb_references = []
    has_plasmabio_paths = False

    # Check each cell
    for i, cell in enumerate(nb.cells):
        if cell.cell_type in ['markdown', 'code']:
            # Check for IFB references
            for pattern in ifb_patterns:
                matches = pattern.findall(cell.source)
                if matches:
                    ifb_references.append(f"Cell {i}: {matches}")

            # Check for Plasmabio paths
            for pattern in plasmabio_patterns:
                if pattern.search(cell.source):
                    has_plasmabio_paths = True

    # Print verification results
    print("\nVerification Results:")
    if ifb_references:
        print("❌ Found IFB references:")
        for ref in ifb_references:
            print(f"  - {ref}")
    else:
        print("✅ No IFB references found")

    print(f"{'✅' if has_plasmabio_paths else '❌'} Plasmabio paths present")

    # Overall verification
    success = (not ifb_references) and has_plasmabio_paths
    print(f"\nOverall verification: {'✅ PASSED' if success else '❌ FAILED'}")
    return success

if __name__ == "__main__":
    verify_plasmabio_compatibility("PS5-2024-merged.ipynb")
