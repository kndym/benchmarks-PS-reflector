"""Compare strict vs relaxed optical stopping rules on equal-count cells.

Runs run_equal_count_cell.py for each (cell, setting) pair in parallel, then
writes a summary of relative cost changes and wall times against the strict
reference (marginal 1e-7, potential change 1e-5, PF tolerance 1e-9).
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
CELLS = [(3600, 1/480), (1600, 1/1280), (1600, 1/320), (400, 1/5120), (400, 1/1280), (400, 1/320), (400, 1/80)]
SETTINGS = [  # (optical marginal L1, potential change, PF EOT tolerance)
    (1e-7, 1e-5, 1e-9), (1e-6, 1e-5, 1e-9), (1e-5, 1e-5, 1e-9), (1e-4, 1e-5, 1e-9),
    (1e-4, 1e-4, 1e-9), (1e-3, 1e-4, 1e-9), (1e-3, 1e-3, 1e-9), (1e-4, 1e-4, 1e-6),
]
METRICS = ["hard_cost", "entropic_cost", "hard_debiased_divergence", "entropic_debiased_divergence"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=5)
    args = parser.parse_args()
    root = args.output_root.resolve()
    root.mkdir(parents=True, exist_ok=True)
    jobs = [(nk, eps, *setting) for nk, eps in CELLS for setting in SETTINGS]

    def run(job):
        nk, eps, marginal, change, pf = job
        tag = f"NK{nk}_eps{eps:.6g}_m{marginal:g}_t{change:g}_pf{pf:g}"
        with (root / f"{tag}.driver.log").open("w") as log:
            code = subprocess.run([sys.executable, "-B", "scripts/run_equal_count_cell.py", "--nk", str(nk),
                                   "--epsilon", repr(eps), "--marginal-tolerance", repr(marginal),
                                   "--refractor-tolerance", repr(change), "--pf-tolerance", repr(pf),
                                   "--output-root", str(root / f"pf{pf:g}")],
                                  cwd=ROOT, stdout=log, stderr=subprocess.STDOUT).returncode
        print(f"{tag}: exit {code}", flush=True)

    with ThreadPoolExecutor(args.jobs) as pool:
        list(pool.map(run, jobs))
    summarize(root)


def summarize(root):
    results = [json.loads(p.read_text()) for p in root.rglob("cell.json")]
    out = []
    for nk, eps in CELLS:
        rows = [r for r in results if r["nk"] == nk and abs(r["epsilon"] - eps) < 1e-15]
        ref = next((r for r in rows if (r["marginal_tolerance"], r["refractor_tolerance"], r["pf_tolerance"])
                    == SETTINGS[0]), None)
        for r in sorted(rows, key=lambda r: (r["pf_tolerance"], r["marginal_tolerance"], r["refractor_tolerance"])):
            entry = dict(nk=nk, epsilon=eps, marginal=r["marginal_tolerance"], change=r["refractor_tolerance"],
                         pf=r["pf_tolerance"], seconds=r["total_seconds"],
                         iterations=[r["final_iterations"], r["source_identity_iterations"],
                                     r["target_identity_iterations"]])
            if ref:
                entry["speedup"] = ref["total_seconds"] / r["total_seconds"]
                entry["relative_change"] = {m: abs(r[m] - ref[m]) / abs(ref[m]) for m in METRICS if m in r}
            out.append(entry)
    (root / "summary.json").write_text(json.dumps(out, indent=2))
    for e in out:
        rel = e.get("relative_change", {})
        print(f"NK={e['nk']:5d} eps={e['epsilon']:.3g} m={e['marginal']:g} t={e['change']:g} pf={e['pf']:g} "
              f"t={e['seconds']/60:6.1f}min x{e.get('speedup', float('nan')):5.2f} "
              + " ".join(f"{k.split('_')[0][0]}{'D' if 'div' in k else 'C'}={v:.2e}" for k, v in rel.items()))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--summarize":
        summarize(Path(sys.argv[2]).resolve())
    else:
        main()
