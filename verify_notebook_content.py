import nbformat
import os
import re

def verify_notebook_content(notebook_path):
    print(f"Verifying notebook: {notebook_path}")

    # Read the notebook
    with open(notebook_path, 'r', encoding='utf-8') as f:
        nb = nbformat.read(f, as_version=4)

    # Initialize verification checklist
    checks = {
        'tool_explanations': False,
        'parameter_docs': False,
        'documentation_links': False,
        'progressive_steps': False,
        'results_organization': False
    }

    # Check markdown cells for explanations and documentation
    markdown_content = ''
    for cell in nb.cells:
        if cell.cell_type == 'markdown':
            content = cell.source.lower()
            markdown_content += content

            # Check for tool explanations
            if any(tool in content for tool in ['qualimap', 'multiqc', 'samtools']):
                checks['tool_explanations'] = True

            # Check for parameter documentation
            if any(param in content for param in ['parameter', 'option', '--']):
                checks['parameter_docs'] = True

            # Check for documentation links
            if 'http' in content or 'www.' in content:
                checks['documentation_links'] = True

            # Check for progressive step explanations
            if any(step in content for step in ['step', 'first', 'then', 'next']):
                checks['progressive_steps'] = True

    # Check code cells for results organization
    for cell in nb.cells:
        if cell.cell_type == 'code':
            if 'Results/' in cell.source and any(tool in cell.source for tool in ['samtools', 'qualimap', 'multiqc']):
                checks['results_organization'] = True

    # Print verification results
    print("\nNotebook Verification Results:")
    print("------------------------------")
    for check, passed in checks.items():
        status = "✓" if passed else "✗"
        print(f"{status} {check.replace('_', ' ').title()}")

    # Additional analysis
    print("\nDetailed Analysis:")
    print("----------------")
    print(f"Total cells: {len(nb.cells)}")
    print(f"Markdown cells: {sum(1 for cell in nb.cells if cell.cell_type == 'markdown')}")
    print(f"Code cells: {sum(1 for cell in nb.cells if cell.cell_type == 'code')}")

    # Suggestions for improvement
    print("\nSuggestions for Improvement:")
    print("--------------------------")
    if not checks['tool_explanations']:
        print("- Add more detailed explanations about tools (qualimap, multiqc, samtools)")
    if not checks['parameter_docs']:
        print("- Include documentation about important parameters and options")
    if not checks['documentation_links']:
        print("- Add links to official documentation and resources")
    if not checks['progressive_steps']:
        print("- Break down complex processes into progressive steps")
    if not checks['results_organization']:
        print("- Ensure results are organized in tool-specific folders")

if __name__ == "__main__":
    notebook_path = "/home/ubuntu/workspace/rnaseq/CEA_rnaseq/temp_repo/PS4-2024-merged.ipynb"
    verify_notebook_content(notebook_path)
