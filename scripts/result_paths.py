"""Repository result locations, with migrated-path and legacy-read compatibility.

Readers accept the old flat bundles. Writers put artifacts in documented,
typed subfolders so a new run cannot scatter files through parent directories.
"""
from datetime import datetime
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
_catalog_stamp = None
_catalog_data = {}


def _catalog():
    global _catalog_stamp, _catalog_data
    path = ROOT / "results/path_catalog.json"
    stamp = path.stat().st_mtime_ns if path.exists() else None
    if stamp != _catalog_stamp:
        _catalog_data = json.loads(path.read_text(encoding="utf-8")) if stamp else {}
        _catalog_stamp = stamp
    return _catalog_data


def resolve_path(path):
    """Resolve known repository aliases, leaving unrelated temporary paths alone."""
    original = Path(str(path).replace("\\", "/"))
    candidate = original if original.is_absolute() else ROOT / original
    try:
        relative = candidate.resolve().relative_to(ROOT.resolve()).as_posix()
    except ValueError:
        return original
    catalog = _catalog()
    exact = catalog.get("files", {}).get(relative)
    if exact:
        return ROOT / exact
    directories = catalog.get("directories", {})
    for old, new in sorted(directories.items(), key=lambda item: len(item[0]), reverse=True):
        old = old.replace("\\", "/").rstrip("/")
        if relative == old or relative.startswith(old + "/"):
            suffix = relative[len(old):].lstrip("/")
            return ROOT / new / suffix
    return original


def ensure_readme(directory, description=None):
    """Create explanatory directory documentation only when it is missing."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    readme = directory / "README.md"
    created_readme = not readme.exists()
    if not readme.exists():
        title = directory.name.replace("_", " ").title()
        descriptions = {
            "data": "Numerical arrays and tabular results (NPZ, NPY, CSV, TSV or TXT). These are the saved values used by the parent run's plots and analysis.",
            "plots": "Rendered figures for this run. Source values live in the parent run's data folder; settings and provenance live in metadata.",
            "metadata": "Experiment parameters, validation reports, provenance and hash manifests. Manifest hash paths are relative to the documented hash_base.",
            "checkpoints": "Reusable numerical checkpoints, separated into optical solutions, push-forward evaluations and self costs.",
            "solutions": "Saved optical solve NPZ checkpoints, including potentials and pushed directions. Reused when a matching run resumes.",
            "pushforward": "Converged push-forward-to-target EOT evaluations and their stopping diagnostics.",
            "self_costs": "Content-hashed self-EOT cache entries used to debias transport costs without repeating optical solves.",
            "diagnostics": "Per-case numerical checks, separated by the quantity being measured.",
            "marginals": "Source and target optical-plan marginal residuals, iteration counts and potential changes.",
            "divergence": "Per-case debiased Sinkhorn evaluations, all three EOT terms, convergence diagnostics and source checkpoint hashes.",
            "logs": "Local process logs, errors and process identifiers. These runtime files are not numerical results.",
            "interactive": "Interactive HTML visualizations generated from saved numerical results.",
        }
        readme.write_text(f"# {title}\n\n" + (description or descriptions.get(directory.name) or
            "Generated numerical artifacts for the parent run. See its README and metadata for settings and provenance.") + "\n",
            encoding="utf-8")
    if created_readme and directory.resolve().is_relative_to(ROOT.resolve()):
        try:
            from check_repository_layout import refresh_directory_readme
        except ModuleNotFoundError:
            from scripts.check_repository_layout import refresh_directory_readme
        refresh_directory_readme(directory)
        if (directory.parent / "README.md").exists():
            refresh_directory_readme(directory.parent)
    return directory


def artifact_group(name):
    name = Path(name).name
    if name.startswith("solve_") and name.endswith(".npz"):
        return Path("checkpoints/solutions")
    if name.startswith("pf_") and name.endswith(".json"):
        return Path("checkpoints/pushforward")
    if name.startswith("residual_") and name.endswith(".json"):
        return Path("diagnostics/marginals")
    if name.startswith("divergence_") and name.endswith(".json"):
        return Path("diagnostics/divergence")
    suffixes = Path(name).suffixes
    suffix = suffixes[-1].lower() if suffixes else ""
    if suffix in (".png", ".jpg", ".jpeg", ".svg", ".pdf"):
        return Path("plots")
    if suffix == ".html":
        return Path("interactive")
    if suffix in (".csv", ".tsv", ".npz", ".npy", ".txt"):
        return Path("data")
    if suffix in (".log", ".err", ".pid"):
        return Path("logs")
    return Path("metadata")


def artifact_path(bundle, name, *, for_write=False):
    """Find a migrated artifact, or return/create its canonical write path."""
    original_bundle = Path(bundle)
    bundle = resolve_path(bundle)
    name = Path(name)
    if name.is_absolute():
        return resolve_path(name)
    legacy = resolve_path(original_bundle / name)
    if len(name.parts) > 1:
        canonical = resolve_path(bundle / name)
    elif bundle.name in {"data", "plots", "metadata", "interactive", "logs",
                         "solutions", "pushforward", "marginals", "divergence", "self_costs"}:
        canonical = bundle / name
    else:
        canonical = bundle / artifact_group(name) / name
    if not for_write:
        for candidate in (canonical, legacy, bundle/name):
            if candidate.exists():
                return candidate
        return canonical
    # A new canonical bundle is documented from its root through every ancestor.
    ensure_readme(bundle)
    try:
        parts = canonical.parent.relative_to(bundle).parts
    except ValueError:
        parts = ()
    parent = bundle
    for part in parts:
        parent = ensure_readme(parent / part)
    return canonical


def artifact_glob(bundle, pattern):
    """Find checkpoint names across typed subfolders and old flat fixtures."""
    bundle = resolve_path(bundle)
    if not bundle.exists():
        return []
    return sorted({resolve_path(p) for p in bundle.rglob(pattern) if p.is_file()})


def dated_run_dir(category, label, *, date=None):
    """Create a descriptive dated run, with documented category parents."""
    date = date or datetime.now().strftime("%b_%d_%y").lower()
    label = re.sub(r"[^a-z0-9_]+", "_", label.lower()).strip("_")
    result_root = ensure_readme(ROOT / "results", "Numerical results organized by experiment family. See path_catalog.json for migrated locations.")
    category_root = result_root
    for part in Path(category).parts:
        category_root = ensure_readme(category_root / part)
    return ensure_readme(category_root / f"{label}_{date}")


def self_cost_dir(bundle, *, for_write=False):
    bundle = resolve_path(bundle)
    canonical = bundle / "checkpoints/self_costs"
    legacy = resolve_path(bundle / "self_cache")
    if not for_write and legacy.exists():
        return legacy
    if for_write:
        ensure_readme(bundle)
        ensure_readme(bundle / "checkpoints")
        ensure_readme(canonical)
    return canonical


def latest_reference_run(pattern, artifact):
    """Find the latest populated reference run without creating a directory."""
    category = ROOT / "results/reference_runs"
    candidates = [(artifact_path(path, artifact), path) for path in category.glob(pattern)
                  if path.is_dir()]
    candidates = [(saved, path) for saved, path in candidates if saved.is_file()]
    return max(candidates, key=lambda item: item[0].stat().st_mtime_ns)[1] if candidates else None
