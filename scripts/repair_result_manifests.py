"""Verify or refresh artifact hashes after a deliberate repository relocation.

Use --apply only after documented path or README changes. Numerical arrays are
never rewritten. Original locations and hashes remain in results/path_catalog.json.
"""
import argparse
import hashlib
import json
import posixpath
from pathlib import Path

try:
    from result_paths import ROOT, artifact_path, resolve_path
except ModuleNotFoundError:
    from scripts.result_paths import ROOT, artifact_path, resolve_path


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()
    catalog = json.loads((ROOT/"results/path_catalog.json").read_text())
    files = catalog["files"]
    reverse = {v:k for k,v in files.items()}
    paths = sorted(p for base in [ROOT/"results", ROOT/"geometric_optics_progress/metadata"]
                   for p in base.rglob("*.json") if p.name != "path_catalog.json")
    unresolved = set()
    verified = 0

    def find(path, value, old_scope=None):
        if not isinstance(value, str):
            return None
        value = value.replace("\\", "/")
        possible = [resolve_path(value), path.parent/value]
        if old_scope:
            old_key = posixpath.normpath(old_scope+"/"+value)
            if old_key in files:
                possible.insert(0, ROOT/files[old_key])
            possible.append(resolve_path(ROOT/old_key))
        for candidate in possible:
            if candidate.exists():
                return candidate.resolve()
        return None

    for iteration in range(8):
        changed = 0
        verified = 0
        for path in paths:
            try:
                original = json.loads(path.read_text(encoding="utf-8-sig"))
            except (json.JSONDecodeError, UnicodeError):
                continue
            obj = json.loads(json.dumps(original))
            old_key = reverse.get(path.relative_to(ROOT).as_posix(), path.relative_to(ROOT).as_posix())
            old_scope = posixpath.dirname(old_key)
            scope = resolve_path(ROOT/old_scope)

            def update(record):
                nonlocal verified
                if isinstance(record, list):
                    for item in record:
                        update(item)
                    return
                if not isinstance(record, dict):
                    return
                for key, value in list(record.items()):
                    if isinstance(value, (dict, list)):
                        update(value)
                for field, sibling in [("cross_checkpoint_sha256", "cross_checkpoint"),
                                       ("saved_map_sha256", "saved_map")]:
                    if field in record and sibling in record:
                        target = find(path, record[sibling], old_scope)
                        if target and target.is_file():
                            record[field] = digest(target); verified += 1
                        else:
                            unresolved.add((str(path.relative_to(ROOT)), str(record[sibling])))
                if "sha256" in record:
                    for sibling in ("archived", "path", "file", "checkpoint"):
                        target = find(path, record.get(sibling), old_scope)
                        if target and target.is_file():
                            record["sha256"] = digest(target); verified += 1
                            break
                for field, candidates in [("npz_sha256", ("npz", "npz_file", "npz_path")),
                                          ("image_sha256", ("image", "image_file", "image_path"))]:
                    if field in record:
                        for sibling in candidates:
                            target = find(path, record.get(sibling), old_scope)
                            if target and target.is_file():
                                record[field] = digest(target); verified += 1
                                break
                for field in ("files_sha256", "artifact_sha256"):
                    mapping = record.get(field)
                    if not isinstance(mapping, dict):
                        continue
                    base_field = "hash_base" if field == "files_sha256" else "artifact_hash_base"
                    existing_base = record.get(base_field)
                    base = (path.parent/existing_base).resolve() if existing_base else scope
                    revised = {}
                    for name, value in mapping.items():
                        target = ((base/name).resolve() if (base/name).is_file() else
                                  find(path, name, old_scope))
                        if target and target.is_file() and target != path:
                            relative = Path(__import__('os').path.relpath(target, base)).as_posix()
                            revised[relative] = digest(target); verified += 1
                        elif target == path:
                            continue  # A manifest cannot contain its own cryptographic hash.
                        else:
                            revised[name] = value
                            unresolved.add((str(path.relative_to(ROOT)), name))
                    record[field] = revised
                    record[base_field] = Path(__import__('os').path.relpath(base, path.parent)).as_posix()
                mapping = record.get("source_checkpoint_sha256")
                source = find(path, record.get("source"), old_scope)
                if isinstance(mapping, dict) and source and source.is_dir():
                    for name in mapping:
                        target = artifact_path(source, name)
                        if target.is_file():
                            mapping[name] = digest(target); verified += 1
                    record["source_checkpoint_hash_base"] = source.relative_to(ROOT).as_posix()
                if isinstance(record.get("csv_sha256"), dict):
                    comparison = find(path, record.get("comparison"), old_scope)
                    if comparison:
                        for mode in record["csv_sha256"]:
                            target = artifact_path(comparison/mode, "costs.csv")
                            if target.is_file():
                                record["csv_sha256"][mode] = digest(target); verified += 1
                if "costs_sha256" in record:
                    target = artifact_path(scope, "costs.csv")
                    if target.is_file():
                        record["costs_sha256"] = digest(target); verified += 1

            update(obj)
            if obj != original:
                changed += 1
                if args.apply:
                    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
        if not changed:
            break
        if not args.apply:
            raise SystemExit(f"{changed} manifests need a documented refresh; run with --apply after intentional relocation.")
    else:
        raise SystemExit("Manifest refresh did not converge; inspect hash dependencies.")
    if unresolved:
        print("Unavailable historical hash targets:")
        for source, name in sorted(unresolved):
            print(f"  {source}: {name}")
        raise SystemExit(1)
    print(f"Verified {verified} artifact hash references across {len(paths)} JSON files.")


if __name__ == "__main__":
    main()
