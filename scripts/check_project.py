#!/usr/bin/env python3
"""Static self-consistency checks for the nc-p3 research project (no network)."""
from pathlib import Path
import json
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
REQUIRED = [
    "README.md", "AGENTS.md", "START_CODEX_ja.md",
    "research/plan_ja.tex", "research/notes_ja.tex",
    "research/BASELINE.md", "research/STATE.json",
    "nc-p3.tex", ".latexmkrc", ".gitignore",
]
errors = []
for rel in REQUIRED:
    if not (ROOT / rel).is_file():
        errors.append("missing: " + rel)

state_path = ROOT / "research/STATE.json"
if state_path.exists():
    state = json.loads(state_path.read_text(encoding="utf-8"))
    if state.get("next_task_id") != "Z4-WINDOW-001":
        errors.append("Unexpected next task; synchronize instructions")
    for path in state.get("document_paths", {}).values():
        if not (ROOT / path).exists():
            errors.append("Broken document path: " + path)
    if len({task["id"] for task in state["tasks"]}) != len(state["tasks"]):
        errors.append("Duplicate task ID")

labels = []
for rel in ["nc-p3.tex", "research/plan_ja.tex", "research/notes_ja.tex"]:
    path = ROOT / rel
    if not path.exists():
        continue
    source = path.read_text(encoding="utf-8")
    depth = 0
    for line in source.splitlines():
        # TeX comment stripping, approximate but catches accidental brace loss
        line = re.split(r"(?<!\\)%", line, maxsplit=1)[0]
        for char in re.findall(r"(?<!\\)[{}]", line):
            depth += 1 if char == "{" else -1
            if depth < 0:
                errors.append("Closing brace without open: " + rel)
                depth = 0
    if depth:
        errors.append("Unbalanced curly braces: " + rel)
    labels += re.findall(r"\\label\{([^}]+)\}", source)
    envs = ["theorem", "lemma", "proposition", "corollary", "proof", "remark", "definition", "document", "abstract"]
    for env in envs:
        if source.count(r"\begin{" + env + "}") != source.count(r"\end{" + env + "}"):
            errors.append("Unbalanced TeX environment " + env + ": " + rel)
if len(labels) != len(set(labels)):
    errors.append("Duplicate TeX labels")

if (ROOT / ".github" / "workflows").exists():
    errors.append("GitHub Actions not allowed")
if errors:
    for err in errors:
        print("FAIL:", err)
    sys.exit(1)
print("PASS: nc-p3 file presence, paths, state, TeX structure (static only)")
print("NOTE: this does not compile TeX or prove any mathematical theorem.")
