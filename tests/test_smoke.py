import json
import re
from pathlib import Path

def test_notebooks_have_quiz_and_discussion():
    # find notebooks recursively under synedu/
    nb_paths = sorted(Path('synedu').rglob('S*.ipynb'))
    assert nb_paths, "No notebooks found under synedu/ (expected pattern synedu/**/S*.ipynb)"

    # regex to match markdown headers like: "# Discussion", "## Discussion", ... (case-insensitive)
    def has_header(text: str, name: str) -> bool:
        pattern = re.compile(r'^\s*#{1,6}\s+' + re.escape(name) + r'\b', re.M | re.I)
        return bool(pattern.search(text))

    for p in nb_paths:
        try:
            nb = json.loads(p.read_text(encoding='utf-8'))
        except Exception as exc:
            raise AssertionError(f"Failed to read/parse notebook {p}: {exc}")

        # concatenate all markdown cells into one text block
        md_cells = [ ''.join(c.get('source', '')) for c in nb.get('cells', []) if c.get('cell_type') == 'markdown' ]
        text = '\n'.join(md_cells)

        if not has_header(text, 'Discussion'):
            raise AssertionError(f"{p} missing a markdown header for 'Discussion' (e.g. '# Discussion' or '## Discussion')")

        if not has_header(text, 'Quiz'):
            raise AssertionError(f"{p} missing a markdown header for 'Quiz' (e.g. '# Quiz' or '## Quiz')")
