import json
import re

def extract_r_dependencies(notebook_path):
    with open(notebook_path, 'r') as f:
        notebook = json.load(f)

    libraries = set()
    programs = set()

    for cell in notebook['cells']:
        if cell['cell_type'] == 'code':
            source = ''.join(cell['source'])

            # Find library() calls
            lib_matches = re.findall(r'library\((.*?)\)', source)
            libraries.update(lib_matches)

            # Find require() calls
            req_matches = re.findall(r'require\((.*?)\)', source)
            libraries.update(req_matches)

            # Find install.packages() calls
            install_matches = re.findall(r'install\.packages\(["\'](.+?)["\']\)', source)
            libraries.update(install_matches)

            # Find BiocManager::install() calls
            bioc_matches = re.findall(r'BiocManager::install\(["\'](.+?)["\']\)', source)
            libraries.update(bioc_matches)

    return {
        'libraries': sorted(list(libraries)),
        'programs': sorted(list(programs))
    }

if __name__ == '__main__':
    notebook_path = 'Pipe_07a-R_intro-to-R.ipynb'
    deps = extract_r_dependencies(notebook_path)

    print(f"\nAnalyzing {notebook_path}...")
    print("\nR Libraries found:")
    for lib in deps['libraries']:
        print(f"- {lib}")

    print("\nPrograms found:")
    for prog in deps['programs']:
        print(f"- {prog}")
