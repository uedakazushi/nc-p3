#!/usr/bin/env python3
"""Offline structural checks, including local TeX inputs; not a proof or compile."""
from pathlib import Path
import json
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
REQUIRED = [
    'README.md', 'AGENTS.md', 'START_CODEX_ja.md', 'nc-p3.tex',
    'sections/elliptic_ore_family.tex', 'sections/ore_complete.tex',
    'sections/isogeny_invariant.tex', 'sections/nonore_family.tex',
    'sections/nonore_period_two.tex',
    'sections/ore_boundary.tex', 'sections/nine_closure.tex', 'sections/four_window.tex',
    'sections/local_moduli.tex', 'sections/algebraic_family.tex', 'sections/frobenius_tensor.tex',
    'computations/frobenius_tensor.py', 'computations/verify_frobenius_tensor.py',
    'computations/verify_quantum_tensor_rank.py', 'computations/verify_symmetric_congruence.py', 'computations/tensor_symmetrization.py', 'computations/verify_rank_two_exclusion.py', 'research/frobenius_tensor_certificate.json',
    'research/classification_branches.json', 'research/plan_ja.tex',
    'research/notes_ja.tex', 'research/BASELINE.md', 'research/STATE.json',
    'research/validation.json', 'research/moduli_certificate.json',
    'research/algebraic_family_certificate.json',
    'research/elliptic_family_certificate.json', 'research/ore_derivation_certificate.json',
    'research/isogeny_certificate.json', 'research/nonore_certificate.json',
    'research/nonore_point_certificate.json',
    'research/ore_boundary_certificate.json', 'research/nine_closure_certificate.json', '.latexmkrc', '.gitignore',
    '.vscode/settings.json', 'requirements.txt', 'computations/pbw_complex.py',
    'computations/verify_moduli_jets.py', 'computations/verify_exterior_trace.py',
    'computations/verify_period_witness.py', 'computations/verify_third_order_samples.py',
    'computations/verify_algebraic_family.py', 'computations/verify_elliptic_family.py', 'computations/verify_ore_derivations.py',
    'computations/verify_isogeny_invariant.py', 'computations/verify_nonore_family.py',
    'computations/verify_nonore_point_geometry.py',
    'computations/verify_ore_boundary.py', 'computations/verify_nine_closure_tangent.py',
]
errors = []
for name in REQUIRED:
    if not (ROOT / name).is_file():
        errors.append('missing ' + name)


def expanded_tex(name, ancestors=()):
    """Expand braced local input/include commands relative to repository root."""
    path = (ROOT / name).resolve()
    if not path.is_relative_to(ROOT) or path in ancestors:
        errors.append('invalid or cyclic TeX input ' + name)
        return ''
    if not path.is_file():
        errors.append('missing TeX input ' + name)
        return ''
    text = path.read_text(encoding='utf-8')
    clean = '\n'.join(re.split(r'(?<!\\)%', line, maxsplit=1)[0]
                      for line in text.splitlines())
    def replace(match):
        child = match.group(1)
        if not child.endswith('.tex'):
            child += '.tex'
        return expanded_tex(child, ancestors + (path,))
    return re.sub(r'\\(?:input|include)\{([^}]+)\}', replace, clean)


try:
    state = json.loads((ROOT / 'research/STATE.json').read_text(encoding='utf-8'))
    for path in state['document_paths'].values():
        if not (ROOT / path).is_file():
            errors.append('broken document path ' + path)
    ids = [item['id'] for item in state['tasks']]
    if len(ids) != len(set(ids)):
        errors.append('duplicate task id')
    if state['next_task_id'] not in ids:
        errors.append('next task not present')
    for path in ROOT.rglob('*.json'):
        if 'build' not in path.parts:
            json.loads(path.read_text(encoding='utf-8'))
except (OSError, ValueError, KeyError) as exc:
    errors.append('JSON/state validation: ' + str(exc))
    state = {}

paper_labels = set()
for name in ['nc-p3.tex', 'research/plan_ja.tex', 'research/notes_ja.tex']:
    clean = expanded_tex(name)
    depth = 0
    for char in re.findall(r'(?<!\\)[{}]', clean):
        depth += 1 if char == '{' else -1
        if depth < 0:
            errors.append('closing brace without open ' + name)
            depth = 0
    if depth:
        errors.append('unbalanced braces ' + name)
    stack = []
    for op, env in re.findall(r'\\(begin|end)\{([^}]+)\}', clean):
        if op == 'begin':
            stack.append(env)
        elif not stack or stack.pop() != env:
            errors.append('environment mismatch ' + name)
    if stack:
        errors.append('unclosed environment ' + name)
    labels = re.findall(r'\\label\{([^}]+)\}', clean)
    if len(labels) != len(set(labels)):
        errors.append('duplicate label ' + name)
    for key in re.findall(r'\\(?:eqref|ref)\{([^}]+)\}', clean):
        if key not in labels:
            errors.append('unresolved ref ' + name + ':' + key)
    bibliography = set(re.findall(r'\\bibitem\{([^}]+)\}', clean))
    for entry in re.findall(r'\\cite(?:\[[^]]*\])?\{([^}]+)\}', clean):
        for key in entry.split(','):
            if key.strip() not in bibliography:
                errors.append('unresolved citation ' + name + ':' + key)
    if name == 'nc-p3.tex':
        paper_labels = set(labels)
for result in state.get('mathematical_results', []):
    proof = result.get('proof', '')
    if proof.startswith('nc-p3.tex:') and proof.split(':', 1)[1] not in paper_labels:
        errors.append('broken proof reference ' + result['id'])
if (ROOT / '.github/workflows').exists():
    errors.append('GitHub Actions forbidden')
if errors:
    print('\n'.join('FAIL: ' + error for error in errors))
    sys.exit(1)
print('PASS: current paths, JSON, tasks, expanded TeX inputs, labels and citations.')
print('NOTE: this is not a compile, formal proof, or higher-order obstruction check.')
