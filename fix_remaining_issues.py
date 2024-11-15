import nbformat
import re

def remove_french_text(text):
    """Replace common French phrases with English equivalents"""
    replacements = {
        r'\ble\b': 'the',
        r'\bla\b': 'the',
        r'\bles\b': 'the',
        r'\bun\b': 'a',
        r'\bune\b': 'a',
        r'\bdes\b': 'the',
        r'\bdu\b': 'of the',
        r'\bau\b': 'to the',
        r'\baux\b': 'to the',
        r'\best\b': 'is',
        r'\bsont\b': 'are',
        r'\bavoir\b': 'have',
        r'\bêtre\b': 'be',
        r'\bdans\b': 'in',
        r'\bpour\b': 'for',
        r'\bavec\b': 'with',
        r'\bsur\b': 'on',
        r'Résultats': 'Results',
        r'Commandes': 'Commands',
        r'Fichiers': 'Files',
        r'Paramètres': 'Parameters'
    }

    for french, english in replacements.items():
        text = re.sub(french, english, text, flags=re.IGNORECASE)
    return text

def remove_ifb_references(text):
    """Remove references to IFB server"""
    ifb_patterns = [
        r'.*\bIFB\b.*\n?',
        r'.*\bifb-.*\n?',
        r'.*core\.cluster\.france-bioinformatique\.fr.*\n?',
        r'.*\bifb core cluster\b.*\n?'
    ]

    for pattern in ifb_patterns:
        text = re.sub(pattern, '', text, flags=re.IGNORECASE)
    return text

def ensure_documentation_links(text):
    """Ensure documentation links are properly formatted"""
    if '](http' in text:  # Already has proper markdown links
        return text

    # Add proper markdown formatting to URLs
    url_pattern = r'(https?://[^\s]+)'
    return re.sub(url_pattern, r'[\1](\1)', text)

def fix_error_handling(code):
    """Ensure proper error handling in code cells"""
    if 'try:' in code and 'except' in code:
        return code

    return f"""try:
    print("Starting process...")  # Progress indicator

{code}

    print("Process completed successfully!")  # Success indicator
except FileNotFoundError as e:
    print(f"Error: Input file not found - {{e}}")
    raise
except Exception as e:
    print(f"Error occurred: {{e}}")
    raise"""

def fix_remaining_issues(notebook_file):
    print(f"Fixing remaining issues in: {notebook_file}")

    # Load notebook
    with open(notebook_file, 'r', encoding='utf-8') as f:
        nb = nbformat.read(f, as_version=4)

    # Process each cell
    for cell in nb.cells:
        if cell.cell_type == 'markdown':
            # Fix French text and IFB references in markdown
            cell.source = remove_french_text(cell.source)
            cell.source = remove_ifb_references(cell.source)
            cell.source = ensure_documentation_links(cell.source)
        elif cell.cell_type == 'code':
            # Fix French comments and add error handling in code
            cell.source = remove_french_text(cell.source)
            if not cell.source.startswith('try:'):
                cell.source = fix_error_handling(cell.source)

    # Save modified notebook
    with open(notebook_file, 'w', encoding='utf-8') as f:
        nbformat.write(nb, f)

    print("Remaining issues fixed successfully")

if __name__ == "__main__":
    fix_remaining_issues("PS5-2024-merged.ipynb")
