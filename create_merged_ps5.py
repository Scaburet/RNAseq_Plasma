import nbformat
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell, new_raw_cell

def read_notebook(path):
    with open(path, 'r', encoding='utf-8') as f:
        return nbformat.read(f, as_version=4)

def create_merged_notebook():
    # Read both notebooks
    pipe06 = read_notebook('/home/ubuntu/attachments/Pipe_06-bash_reads-counts-pseudomapping.ipynb')
    ps5 = read_notebook('/home/ubuntu/attachments/PS5-ReadCounts-bash-2022-executed.ipynb')

    # Create new notebook
    merged = new_notebook()

    # Initialize cell counter
    code_cell_counter = 1

    # Add cells from Pipe_06, maintaining its layout
    for cell in pipe06.cells:
        if cell.cell_type == 'code':
            # Update code cell numbering
            source = cell.source
            if not source.startswith('## Code cell'):
                source = f'## Code cell {code_cell_counter} ##\n\n{source}'
            else:
                source = f'## Code cell {code_cell_counter} ##\n\n' + '\n'.join(source.split('\n')[1:])
            new_cell = new_code_cell(source=source)
            new_cell.metadata = cell.metadata
            merged.cells.append(new_cell)
            code_cell_counter += 1
        elif cell.cell_type == 'raw':
            new_cell = new_raw_cell(source=cell.source)
            new_cell.metadata = cell.metadata
            merged.cells.append(new_cell)
        else:  # markdown
            new_cell = new_markdown_cell(source=cell.source)
            new_cell.metadata = cell.metadata
            merged.cells.append(new_cell)

    # Write the merged notebook
    output_path = '/home/ubuntu/workspace/rnaseq/RNAseq_Plasma/PS5-2024-merged.ipynb'
    with open(output_path, 'w', encoding='utf-8') as f:
        nbformat.write(merged, f)

if __name__ == '__main__':
    create_merged_notebook()
