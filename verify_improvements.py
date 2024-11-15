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

    # Check each cell
    for cell in nb.cells:
        # Check markdown cells for documentation features
        if cell.cell_type == 'markdown':
            content = cell.source.lower()

            if 'learning objective' in content or 'objectives:' in content:
                requirements['learning_objectives'] = True

            if 'documentation:' in content or 'reference:' in content:
                requirements['documentation_links'] = True

            if 'parameter' in content and ':' in content:
                requirements['parameter_descriptions'] = True

            if 'troubleshoot' in content or 'common error' in content:
                requirements['troubleshooting'] = True

            if any(re.search(pattern, content) for pattern in french_patterns):
                requirements['english_text'] = False

            if re.match(r'^#+\s+\w+', cell.source):  # Proper header format
                requirements['clear_headers'] = True

        # Check code cells for technical features
        elif cell.cell_type == 'code':
            if 'try:' in cell.source and 'except:' in cell.source:
                requirements['error_handling'] = True

            if 'validate' in cell.source.lower() or 'check' in cell.source.lower():
                requirements['validation_steps'] = True

            if 'progress' in cell.source.lower() or 'status' in cell.source.lower():
                requirements['progress_indicators'] = True

            if 'memory' in cell.source.lower() or 'cpu' in cell.source.lower():
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
