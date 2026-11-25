"""Release assembly with content checks (manual Section 15; acceptance J2-J3).

Refuses to assemble an archive if any project file contains placeholders, tool attribution,
time-based planning language, delegation language, or emoji. Writes checksums for all archive
members and a reproducibility manifest.
"""

from __future__ import annotations

import hashlib
import json
import platform
import re
import subprocess
import sys
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TEXT_EXT = {".py", ".md", ".toml", ".yaml", ".yml", ".cff", ".json", ".txt", ".tex", ".bib", ".cfg"}
EXCLUDE_DIRS = {".git", "data/raw", "dist", "build", ".pytest_cache", ".ruff_cache", ".mypy_cache", "__pycache__", ".hypothesis", ".venv"}
EXCLUDE_FILES = {"make_release.py", "CLAUDE.md"}

PATTERNS = {
    "placeholder": re.compile(r"\b(TODO|TBD|FIXME|XXX|lorem ipsum)\b|\[insert|\[placeholder", re.I),
    "delegation": re.compile(r"\b(you should now|the user will|the researcher should|manually verify|run this yourself|ask the researcher|provide the missing dataset)\b", re.I),
    "schedule": re.compile(r"\b(week \d|weeks? \d|month \d|sprint|timeline|deadlines?(?!\s*=)|milestone \d|day \d+ plan|session notes|dev diary)\b", re.I),
    "tool_attribution": re.compile(r"\b(assisted by|built with|generated with|co-authored-by)\b", re.I),
    "emoji": re.compile("[\U0001F300-\U0001FAFF☀-➿]"),
}


def iter_files() -> list[Path]:
    out = []
    for p in ROOT.rglob("*"):
        if not p.is_file():
            continue
        rel = p.relative_to(ROOT)
        if any(part in EXCLUDE_DIRS or str(rel).startswith(d) for d in EXCLUDE_DIRS for part in rel.parts):
            continue
        if p.suffix.lower() in TEXT_EXT:
            out.append(p)
    return out


def content_check(files: list[Path]) -> list[str]:
    problems = []
    for f in files:
        if f.name in EXCLUDE_FILES:
            continue
        text = f.read_text(errors="ignore")
        for name, pat in PATTERNS.items():
            m = pat.search(text)
            if m:
                problems.append(f"{f.relative_to(ROOT)}: {name}: {m.group(0)!r}")
    return problems


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    h.update(p.read_bytes())
    return h.hexdigest()


def main(version: str) -> int:
    files = iter_files()
    problems = content_check(files)
    if problems:
        print("content check failed:")
        for p in problems:
            print("  ", p)
        return 1
    rel_dir = ROOT / "release"
    rel_dir.mkdir(exist_ok=True)
    archive = rel_dir / f"sitewindow-{version}.tar.gz"
    members = [f for f in files if f.name not in EXCLUDE_FILES] + [ROOT / "LICENSE"]
    with tarfile.open(archive, "w:gz") as tar:
        for f in sorted(set(members)):
            tar.add(f, arcname=str(Path(f"sitewindow-{version}") / f.relative_to(ROOT)))
    checksums = {str(f.relative_to(ROOT)): sha256(f) for f in sorted(set(members))}
    (rel_dir / "CHECKSUMS.sha256").write_text("\n".join(f"{v}  {k}" for k, v in checksums.items()) + "\n")
    manifest = {
        "version": version,
        "archive": archive.name,
        "archive_sha256": sha256(archive),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "git_describe": subprocess.run(["git", "describe", "--always", "--dirty"], capture_output=True, text=True, cwd=ROOT).stdout.strip(),
        "sources_lock": json.loads((ROOT / "data" / "sources.lock").read_text()) if (ROOT / "data" / "sources.lock").exists() else None,
    }
    (rel_dir / "reproducibility_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"assembled {archive} ({len(members)} files); checksums and manifest written")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "0.1.0"))
