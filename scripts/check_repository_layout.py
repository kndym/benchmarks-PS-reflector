"""Refresh directory indexes and check the repository's output organization.

Run with --refresh-readmes after adding or relocating result folders. Source
directories may contain many scripts; result parents must have at most five
regular files. The numerical library and tests are deliberately excluded.
"""
import argparse
import json
import os
from pathlib import Path
import re
from urllib.parse import quote, unquote

ROOT = Path(__file__).resolve().parents[1]
EXCLUDED = {".git", ".research", ".claude", ".cursor", ".agents", ".codex",
            "__pycache__", ".pytest_cache", "refracter", "tests"}
START = "<!-- directory-index:start -->"
END = "<!-- directory-index:end -->"
PURPOSES = {
    "data": "Numerical arrays and tabular results for the parent run. NPZ/NPY files retain their original array names; CSV files retain their numerical columns.",
    "plots": "Static figures for the parent run or archived export. See the parent README and parameter metadata for the metric, units, and experiment settings.",
    "interactive": "Interactive HTML visualizations. Open the HTML file in a browser; local assets are linked relative to this directory.",
    "metadata": "Run settings, validation reports, and provenance records. Hash manifests specify their relative hash base where needed.",
    "checkpoints": "Saved numerical state, separated by solve type. These files allow evaluation or resumption without recomputing completed optical solves.",
    "solutions": "Optical solution NPZ checkpoints: source/target arrays, potentials, and pushed directions. Original numerical filenames and bytes are preserved.",
    "pushforward": "Independent push-forward-to-target entropic transport evaluations, including the objective and recorded convergence diagnostics.",
    "self_costs": "Content-addressed self-transport evaluations reused when forming debiased Sinkhorn divergence. These are caches, not additional optical runs.",
    "diagnostics": "Numerical diagnostics, separated from plotted transport costs and optical solution arrays.",
    "marginals": "Source/target marginal residuals and optical stopping diagnostics. Potential-change stopping and marginal feasibility are distinct checks.",
    "divergence": "Per-case debiased Sinkhorn records, including cross/self terms, convergence diagnostics, source checkpoint paths, and hashes.",
    "logs": "Local process logs and runtime identifiers. These are operational records and are excluded from GitHub.",
    "helpers": "Archived local supervision helpers. These are implementation records, not numerical results or current entry points.",
    "source": "Supporting source documents and historical notes retained for provenance.",
    "packages": "Archived export packages. Canonical editable sources and current numerical data are linked from the parent README.",
}


def directories():
    for base, children, _ in os.walk(ROOT):
        children[:] = sorted(n for n in children if n not in EXCLUDED and not n.startswith('.'))
        yield Path(base)


def description(path):
    if path.name in PURPOSES:
        return PURPOSES[path.name]
    relative = path.relative_to(ROOT).as_posix()
    if "case_" in path.name:
        return "One representative refractor case. Its data, figures, and shared settings are separated below."
    if "smoke_checks" in relative:
        return "Archived pipeline smoke-check output. These small checks are not production numerical experiments."
    if "progress_captures" in relative:
        return "A preserved intermediate progress capture. Its metadata records the actual completed-case counts; newer complete results are under results/published."
    if path.name.endswith("_26"):
        return "A named run, capture, or archived update. The date suffix uses the recorded run/publication date when known; otherwise it records the October 27 catalog update. See metadata and the parent index for provenance."
    if relative.startswith("BenchmarkCode"):
        return "C++ reference implementation sources and fixed sampling inputs. Generated numerical outputs belong under results/reference_runs, not beside these source files."
    if relative.startswith("scripts"):
        return "Python and PowerShell entry points for generation, plotting, validation, and maintenance. Shared output routing is implemented in result_paths.py."
    return "Supporting files grouped by purpose. Follow the directory links below and the parent README for context."


def refresh_directory_readme(path):
    """Refresh one generated index while preserving the explanatory prose."""
    path = Path(path)
    if path == ROOT or not path.is_dir():
        return
    readme = path / "README.md"
    previous = readme.read_text(encoding="utf-8") if readme.exists() else ""
    existing = previous or ("# " + path.name.replace("_", " ").title()
                            + "\n\n" + description(path) + "\n")
    existing = re.sub(re.escape(START)+r".*?"+re.escape(END), "", existing, flags=re.S).rstrip()
    children = sorted(p for p in path.iterdir() if p.is_dir() and p.name not in EXCLUDED)
    lines = [START, "", "## Directory guide", "", "[Parent directory](../README.md)"]
    if children:
        lines += ["", "| Folder | Contents |", "| --- | --- |"]
        lines += [f"| [{p.name}/]({quote(p.name)}/README.md) | {description(p)} |" for p in children]
    elif path.name in PURPOSES:
        lines += ["", PURPOSES[path.name]]
    lines += ["", END]
    revised = existing+"\n\n"+"\n".join(lines)+"\n"
    if revised != previous:
        readme.write_text(revised, encoding="utf-8")


def refresh_readmes():
    for path in directories():
        refresh_directory_readme(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh-readmes", action="store_true")
    args = parser.parse_args()
    if args.refresh_readmes:
        refresh_readmes()
    errors = []
    checked = 0
    for path in directories():
        checked += 1
        if not (path/"README.md").exists():
            errors.append(f"Missing README: {path.relative_to(ROOT)}")
        children = [p for p in path.iterdir() if p.is_dir() and p.name not in EXCLUDED]
        files = [p for p in path.iterdir() if p.is_file()]
        # Script/source roots are the user's stated exception; numerical parents are not.
        source_exception = path == ROOT/"scripts" or path == ROOT/"scripts/archived"
        if children and len(files) > 5 and not source_exception:
            errors.append(f"Crowded parent: {path.relative_to(ROOT)} ({len(files)} files)")
        for md in path.glob("*.md"):
            for target in re.findall(r"!?\[[^\]\n]*\]\(([^)\n]+)\)", md.read_text(encoding="utf-8")):
                target = target.strip().strip('<>')
                if re.match(r"^[a-zA-Z][a-zA-Z0-9+.-]*:", target) or target.startswith('#'):
                    continue
                target = unquote(target.split('#',1)[0])
                if target and not (md.parent/target).exists():
                    errors.append(f"Broken link: {md.relative_to(ROOT)} -> {target}")
    if errors:
        print("\n".join(errors))
        raise SystemExit(1)
    print(f"Layout checked: {checked} directories, README coverage, parent file limits, and local Markdown links passed.")


if __name__ == "__main__":
    main()
