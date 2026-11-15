#!/usr/bin/env python3
"""Offline structural checks; not a mathematical proof or a TeX compile."""
from pathlib import Path
import json
import re
import sys

ROOT=Path(__file__).resolve().parents[1]
required=['README.md','AGENTS.md','START_CODEX_ja.md','nc-p3.tex',
          'research/plan_ja.tex','research/notes_ja.tex','research/BASELINE.md',
          'research/STATE.json','research/validation.json','research/moduli_certificate.json',
          '.latexmkrc','.gitignore','.vscode/settings.json','requirements.txt',
          'computations/pbw_complex.py','computations/verify_moduli_jets.py',
          'computations/verify_exterior_trace.py','computations/verify_period_witness.py',
          'computations/verify_third_order_samples.py',
          'computations/verify_algebraic_family.py','research/algebraic_family_certificate.json']
errors=[]
for name in required:
    if not (ROOT/name).is_file():errors.append('missing '+name)
state=json.loads((ROOT/'research/STATE.json').read_text())
for path in state['document_paths'].values():
    if not (ROOT/path).is_file():errors.append('broken document path '+path)
ids=[x['id'] for x in state['tasks']]
if len(ids)!=len(set(ids)):errors.append('duplicate task id')
if state['next_task_id'] not in ids:errors.append('next task not present')
for path in ROOT.rglob('*.json'):
    if 'build' not in path.parts:json.loads(path.read_text())
labels=set()
for name in ['nc-p3.tex','research/plan_ja.tex','research/notes_ja.tex']:
    source=(ROOT/name).read_text()
    clean='\n'.join(re.split(r'(?<!\\)%',line,maxsplit=1)[0] for line in source.splitlines())
    depth=0
    for ch in re.findall(r'(?<!\\)[{}]',clean):
        depth+=1 if ch=='{' else -1
        if depth<0:errors.append('closing brace without open '+name);depth=0
    if depth:errors.append('unbalanced braces '+name)
    stack=[]
    for op,env in re.findall(r'\\(begin|end)\{([^}]+)\}',clean):
        if op=='begin':stack.append(env)
        elif not stack or stack.pop()!=env:errors.append('environment mismatch '+name)
    if stack:errors.append('unclosed environment '+name)
    labs=re.findall(r'\\label\{([^}]+)\}',clean)
    if len(labs)!=len(set(labs)):errors.append('duplicate label '+name)
    for key in re.findall(r'\\(?:eqref|ref)\{([^}]+)\}',clean):
        if key not in labs:errors.append('unresolved ref '+name+':'+key)
    bib=set(re.findall(r'\\bibitem\{([^}]+)\}',clean))
    for entry in re.findall(r'\\cite(?:\[[^]]*\])?\{([^}]+)\}',clean):
        for key in entry.split(','):
            if key.strip() not in bib:errors.append('unresolved citation '+name+':'+key)
    if name=='nc-p3.tex':labels=set(labs)
for result in state['mathematical_results']:
    if result.get('proof','').startswith('nc-p3.tex:'):
        if result['proof'].split(':',1)[1] not in labels:errors.append('broken proof reference '+result['id'])
if (ROOT/'.github/workflows').exists():errors.append('GitHub Actions forbidden')
if errors:
    print('\n'.join('FAIL: '+x for x in errors));sys.exit(1)
print('PASS: current paths, JSON, task ids, proof labels, citations, TeX braces/environments.')
print('NOTE: this is not a compile, formal proof, or verification of higher-order obstructions.')
