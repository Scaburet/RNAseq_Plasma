import nbformat
import re

def remove_ifb_references(filename):
    print(f"Removing IFB references from: {filename}")

    # Load notebook
    with open(filename, 'r', encoding='utf-8') as f:
        nb = nbformat.read(f, as_version=4)

    # Patterns to identify IFB references
    ifb_patterns = [
        r'.*\bIFB\b.*',
        r'.*\bifb\b.*',
        r'.*\bifb-core\b.*',
        r'.*\bifb\..*',
        r'/shared/.*',
        r'/ifb/.*',
        r'.*\.ifb\..*'
    ]

    # Compile patterns
    patterns = [re.compile(pattern, re.IGNORECASE) for pattern in ifb_patterns]

    # Process each cell
    cells_to_remove = []
    for i, cell in enumerate(nb.cells):
        if cell.cell_type in ['markdown', 'code']:
            # Check if cell contains IFB references
            has_ifb = any(pattern.search(cell.source) for pattern in patterns)

            if has_ifb:
                print(f"Found IFB reference in cell {i}")
                if cell.cell_type == 'markdown':
                    # For markdown, remove only the lines with IFB references
                    new_lines = []
                    for line in cell.source.split('\n'):
                        if not any(pattern.search(line) for pattern in patterns):
                            new_lines.append(line)
                    cell.source = '\n'.join(new_lines)
                else:
                    # For code cells, check if it's entirely IFB-specific
                    if all(pattern.search(cell.source) for pattern in patterns):
                        cells_to_remove.append(i)
                    else:
                        # Remove only IFB-specific lines
                        new_lines = []
                        for line in cell.source.split('\n'):
                            if not any(pattern.search(line) for pattern in patterns):
                                new_lines.append(line)
                        cell.source = '\n'.join(new_lines)

    # Remove cells marked for removal (in reverse order to maintain indices)
    for i in reversed(cells_to_remove):
        del nb.cells[i]

    # Save modified notebook
    with open(filename, 'w', encoding='utf-8') as f:
        nbformat.write(nb, f)

    print("IFB references removal complete")

if __name__ == "__main__":
    remove_ifb_references("PS5-2024-merged.ipynb")
