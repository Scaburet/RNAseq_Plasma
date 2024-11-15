import nbformat

def add_error_handling(code):
    """Add proper error handling to code cells"""
    # Skip cells that already have error handling
    if 'try:' in code and 'except' in code:
        return code

    # Skip cells that are just comments or empty
    if not code.strip() or code.strip().startswith('#'):
        return code

    # Skip raw cells (indicated by comments about being silent)
    if '## Raw Cell' in code or 'silent raw cell' in code.lower():
        return code

    return f"""try:
    # Progress indicator
    print("Starting process...")

    {code.strip()}

    # Success indicator
    print("Process completed successfully!")
except FileNotFoundError as e:
    print(f"Error: Input file not found - {{e}}")
    raise
except Exception as e:
    print(f"Error occurred: {{e}}")
    raise"""

def fix_notebook_error_handling(notebook_file):
    print(f"Adding error handling to: {notebook_file}")

    # Load notebook
    with open(notebook_file, 'r', encoding='utf-8') as f:
        nb = nbformat.read(f, as_version=4)

    # Process each cell
    for cell in nb.cells:
        if cell.cell_type == 'code':
            cell.source = add_error_handling(cell.source)

    # Save modified notebook
    with open(notebook_file, 'w', encoding='utf-8') as f:
        nbformat.write(nb, f)

    print("Error handling added successfully")

if __name__ == "__main__":
    fix_notebook_error_handling("PS5-2024-merged.ipynb")
