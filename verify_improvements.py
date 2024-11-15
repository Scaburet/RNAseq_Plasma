import nbformat
import re

def verify_improvements(notebook_file):
    print(f"Verifying improvements in: {notebook_file}")

    # Load notebook
    with open(notebook_file, 'r', encoding='utf-8') as f:
        nb = nbformat.read(f, as_version=4)

    # Required sections and features to check
    requirements = {
        'learning_objectives': False,
        'documentation_links': False,
        'parameter_descriptions': False,
        'error_handling': False,
        'validation_steps': False,
        'troubleshooting': False,
        'english_text': True,  # Assume English until French is found
        'clear_headers': False,
        'progress_indicators': False,
        'resource_requirements': False
    }

    # French text patterns
    french_patterns = [
        r'\b(le|la|les|un|une|des|du|au|aux)\b',
        r'\b(est|sont|avoir|être)\b',
        r'\b(dans|pour|avec|sur)\b'
    ]

    # Documentation link pattern
    doc_link_pattern = r'\[([^\]]+)\]\((https?://[^\s)]+)\)'

    # Check each cell
    for cell in nb.cells:
        cell_content = cell.source.lower()

        # Check markdown cells for documentation features
        if cell.cell_type == 'markdown':
            if 'learning objective' in cell_content or 'objectives:' in cell_content:
                requirements['learning_objectives'] = True

            if re.search(doc_link_pattern, cell.source, re.IGNORECASE):
                requirements['documentation_links'] = True

            if 'parameter' in cell_content and ':' in cell_content:
                requirements['parameter_descriptions'] = True

            if 'troubleshoot' in cell_content or 'common error' in cell_content:
                requirements['troubleshooting'] = True

            if any(re.search(pattern, cell_content) for pattern in french_patterns):
                requirements['english_text'] = False

            if re.match(r'^#+\s+\w+', cell.source):  # Proper header format
                requirements['clear_headers'] = True

        # Check code cells for technical features
        elif cell.cell_type == 'code':
            # Check for error handling with specific exceptions
            if ('try:' in cell.source and
                ('except FileNotFoundError' in cell.source or
                 'except Exception' in cell.source or
                 'except:' in cell.source)):
                requirements['error_handling'] = True

            if 'validate' in cell_content or 'check' in cell_content:
                requirements['validation_steps'] = True

            if 'print("starting' in cell_content or 'print("process' in cell_content:
                requirements['progress_indicators'] = True

            if 'memory' in cell_content or 'cpu' in cell_content or 'thread' in cell_content:
                requirements['resource_requirements'] = True

    # Print verification results
    print("\nVerification Results:")
    for req, status in requirements.items():
        print(f"{'✅' if status else '❌'} {req.replace('_', ' ').title()}")

    # Overall verification
    success = all(requirements.values())
    print(f"\nOverall verification: {'✅ PASSED' if success else '❌ FAILED'}")
    return success

if __name__ == "__main__":
    verify_improvements("PS5-2024-merged.ipynb")
