import nbformat
from typing import Dict, List, Tuple

def count_cells_by_type(notebook_path: str) -> Dict[str, int]:
    """Count the number of cells by type in the notebook."""
    with open(notebook_path, 'r', encoding='utf-8') as f:
        nb = nbformat.read(f, as_version=4)

    counts = {'markdown': 0, 'code': 0, 'raw': 0}
    for cell in nb.cells:
        counts[cell.cell_type] = counts.get(cell.cell_type, 0) + 1
    return counts

def analyze_markdown_content(notebook_path: str) -> Dict[str, bool]:
    """Analyze markdown cells for required content."""
    with open(notebook_path, 'r', encoding='utf-8') as f:
        nb = nbformat.read(f, as_version=4)

    analysis = {
        'has_step_explanations': False,
        'has_file_requirements': False,
        'has_parameter_docs': False,
        'has_documentation_links': False,
        'has_genome_indexing': False
    }

    for cell in nb.cells:
        content = cell.source.lower()

        # Check markdown cells
        if cell.cell_type == 'markdown':
            if any(word in content for word in ['step', 'phase', 'stage']):
                analysis['has_step_explanations'] = True
            if any(phrase in content for phrase in ['required file', 'input file', 'required input']):
                analysis['has_file_requirements'] = True
            if any(word in content for word in ['parameter', 'option', '--']):
                analysis['has_parameter_docs'] = True
            if 'http' in content or 'documentation' in content:
                analysis['has_documentation_links'] = True
            if 'genome' in content and ('index' in content or 'reference' in content):
                analysis['has_genome_indexing'] = True

        # Check raw cells for genome indexing commands
        elif cell.cell_type == 'raw':
            if ('wget' in content and 'genome' in content) or ('samtools faidx' in content):
                analysis['has_genome_indexing'] = True

        # Check code cells for genome indexing commands
        elif cell.cell_type == 'code':
            if ('wget' in content and 'genome' in content) or ('samtools faidx' in content):
                analysis['has_genome_indexing'] = True

    return analysis

def check_results_structure(notebook_path: str) -> Dict[str, bool]:
    """Check if the notebook creates proper Results folder structure."""
    with open(notebook_path, 'r', encoding='utf-8') as f:
        nb = nbformat.read(f, as_version=4)

    analysis = {
        'has_results_folder': False,
        'has_tool_subfolders': False,
        'processes_first_two_samples': False
    }

    for cell in nb.cells:
        if cell.cell_type in ['code', 'raw']:
            content = cell.source.lower()
            if 'mkdir' in content and 'results' in content:
                analysis['has_results_folder'] = True
            if any(f'results/{tool}' in content.replace(' ', '') for tool in ['samtools', 'qualimap', 'multiqc']):
                analysis['has_tool_subfolders'] = True
            if ('head -n 2' in content or 'first two samples' in content) and ('fastq' in content or '.bam' in content):
                analysis['processes_first_two_samples'] = True

    return analysis

def main():
    notebook_path = 'PS4-2024-merged.ipynb'

    print("Analyzing notebook content...")
    print("\nCell counts:")
    counts = count_cells_by_type(notebook_path)
    for cell_type, count in counts.items():
        print(f"- {cell_type}: {count}")

    print("\nContent analysis:")
    content_analysis = analyze_markdown_content(notebook_path)
    for feature, present in content_analysis.items():
        print(f"- {feature.replace('_', ' ').title()}: {'✓' if present else '✗'}")

    print("\nResults structure:")
    structure_analysis = check_results_structure(notebook_path)
    for feature, present in structure_analysis.items():
        print(f"- {feature.replace('_', ' ').title()}: {'✓' if present else '✗'}")

if __name__ == '__main__':
    main()
