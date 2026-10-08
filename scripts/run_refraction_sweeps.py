"""Four parameter sweeps using the original uniform refractor pipeline.

Run from the repository root: python scripts/run_refraction_sweeps.py
"""
import argparse
import csv
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from scipy.special import logsumexp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from generate_results import (
    REPO_ROOT, SRC_THETA, SRC_PHI, TGT_THETA, TGT_PHI,
    REGULARIZATION_MULTIPLIER, gen_spherical_patch, build_refraction_artifacts,
    P_refraction_patch, Q_refraction_patch,
)
from refracter.cost import set_kappa, cost_matrix_chunk
from refracter.sinkhorn import run_sinkhorn_divergence, entropic_ot_cost
import refracter.sinkhorn as sinkhorn_module


def joint_epsilon(nk):
    return 1.0 / (REGULARIZATION_MULTIPLIER * int(np.floor(np.sqrt(nk))))


def enable_threaded_reductions(workers):
    """Split independent outputs, still evaluating each with SciPy logsumexp."""
    if workers <= 1:
        return
    pool = ThreadPoolExecutor(max_workers=workers)
    original = sinkhorn_module.logsumexp

    def threaded(a, axis=None, **kwargs):
        if a.ndim != 2 or axis not in (0, 1) or a.size < 1_000_000 or kwargs:
            return original(a, axis=axis, **kwargs)
        # Columns are independent for axis=0; rows are independent for axis=1.
        pieces = np.array_split(a, min(workers, a.shape[1-axis]), axis=1-axis)
        return np.concatenate(list(pool.map(lambda piece: original(piece, axis=axis), pieces)))

    sinkhorn_module.logsumexp = threaded


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    grids = parser.add_mutually_exclusive_group()
    grids.add_argument("--expanded", action="store_true", help="20 logarithmic values per axis, spanning 100x")
    grids.add_argument("--extended", action="store_true", help="40 values: NK 100-20000, epsilon 0.0000625-0.025")
    parser.add_argument("--experiments", nargs="+", type=int, choices=[1, 2, 3, 4], default=[1, 2, 3, 4])
    parser.add_argument("--import-costs", type=Path, help="CSV from another machine for experiments 1-3")
    parser.add_argument("--skip-pf-nk", nargs="*", type=int, default=[], help="PF cases assigned to another machine")
    parser.add_argument("--nk", nargs="+", type=int)
    parser.add_argument("--epsilon", nargs="+", type=float)
    parser.add_argument("--pf-epsilon", nargs="+", type=float)
    parser.add_argument("--fixed-nk", type=int, default=1600)
    parser.add_argument("--target-nk", type=int, default=1600, help="Target point count (default 1600); NK varies source only")
    parser.add_argument("--fixed-epsilon", type=float, default=1/320)
    parser.add_argument("--cost-epsilon", type=float, default=1/320,
                        help="Fixed PF-to-target EOT evaluation epsilon for experiments 1-3")
    parser.add_argument("--joint-only", action="store_true",
                        help="Experiment 4 uses only the coupled refractor epsilon(NK) schedule")
    parser.add_argument("--pf-max-iter", type=int, default=2000)
    parser.add_argument("--refractor-max-iter", type=int, default=17, help="Final transport updates; 0 disables cap")
    parser.add_argument("--identity-max-iter", type=int, default=17, help="Final updates per self solve; 0 disables cap")
    parser.add_argument("--refractor-tolerance", type=float, default=1e-5)
    parser.add_argument("--pf-tolerance", type=float, default=1e-9)
    parser.add_argument("--chunk-size", type=int, default=512)
    parser.add_argument("--pf-cache-entries", type=int, default=4_000_000)
    parser.add_argument("--workers", type=int, default=1, help="Threads for independent SciPy reductions")
    parser.add_argument("--pf-warm-start", action="store_true",
                        help="Reuse PF potentials from the preceding larger epsilon")
    parser.add_argument("--output", type=Path, default=Path(REPO_ROOT) / "results/refraction_sweeps")
    args = parser.parse_args()
    if min(args.refractor_max_iter, args.identity_max_iter) < 0:
        parser.error("Iteration caps must be nonnegative (0 means uncapped)")
    if args.chunk_size < 1 or not all(np.isfinite(v) and v > 0 for v in (args.refractor_tolerance, args.pf_tolerance)):
        parser.error("Chunk size and tolerances must be finite and positive")
    enable_threaded_reductions(args.workers)
    default_nk = np.rint(np.geomspace(100, 10000, 20)).astype(int).tolist() if args.expanded else [100, 225, 400, 900, 1600]
    default_epsilon = np.geomspace(.000125, .0125, 20).tolist() if args.expanded else [1/80, 1/160, 1/320, 1/640]
    if args.extended:
        default_nk = np.rint(np.geomspace(100, 20000, 40)).astype(int).tolist()
        default_epsilon = np.geomspace(.0000625, .025, 40).tolist()
    args.nk = args.nk or default_nk
    args.epsilon = args.epsilon or default_epsilon
    args.pf_epsilon = args.pf_epsilon or default_epsilon
    args.nk = sorted(set(args.nk))
    if min(args.nk + [args.fixed_nk]) < 2:
        parser.error("NK must be at least 2")
    if args.target_nk is not None and args.target_nk < 2:
        parser.error("target NK must be at least 2")
    if args.target_nk is not None and args.import_costs:
        parser.error("Fixed-target runs cannot import legacy costs without target-size metadata")
    if any(not 0 < e <= 1 for e in args.epsilon + args.pf_epsilon + [args.fixed_epsilon, args.cost_epsilon]):
        parser.error("epsilon values must be in (0, 1]")
    args.output.mkdir(parents=True, exist_ok=True)
    args.self_cost_mode = "separate_marginals_v1" if args.target_nk is not None else "legacy"
    args.cost_metric = "pushforward_target_eot_squared_euclidean"
    (args.output / "parameters.json").write_text(json.dumps(vars(args), default=str, indent=2), encoding="utf-8")
    set_kappa(0.6)
    cache = {}
    solver_tag = ""
    if (args.refractor_max_iter, args.identity_max_iter, args.refractor_tolerance) != (17, 17, 1e-5):
        solver_tag = f"_iter{args.refractor_max_iter}_id{args.identity_max_iter}_tol{args.refractor_tolerance:.17g}"
    rows = []
    imported = {}
    imported_pf = {}
    if args.import_costs:
        with args.import_costs.open(newline="", encoding="utf-8") as stream:
            for row in csv.DictReader(stream):
                key = (int(row["experiment"]), int(row["nk"]), float(row["refractor_epsilon"]))
                if row["experiment"] == "4":
                    imported_pf[key[1], key[2], float(row["pf_epsilon"])] = row
                else:
                    imported[key] = row

    def write_costs():
        with (args.output / "costs.csv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)

    def solve(nk, epsilon):
        key = (nk, epsilon)
        target_nk = args.target_nk or nk
        target_tag = f"_targetNK{target_nk}_self" if args.target_nk is not None else ""
        target_tag += solver_tag
        if key not in cache:
            checkpoint = args.output / f"solve_NK{nk}{target_tag}_eps{epsilon:.17g}.npz"
            if checkpoint.exists():
                with np.load(checkpoint) as data:
                    cache[key] = ({"total_cost": float(data["total_cost"])},
                                  {"y_push_3d": data["y_push_3d"]},
                                  data["y"], data["p"], data["q"])
                return cache[key]
            print(f"Refractor NK={nk}, epsilon={epsilon:.8g}", flush=True)
            x = gen_spherical_patch(nk, *SRC_THETA, *SRC_PHI, skip=0)
            y = gen_spherical_patch(target_nk, *TGT_THETA, *TGT_PHI, skip=0)
            p, q = P_refraction_patch(x), Q_refraction_patch(y)
            p, q = p/p.sum(), q/q.sum()
            solution = run_sinkhorn_divergence(x, y, p, q, verbose=False, epsilon=epsilon,
                                               separate_marginals=args.target_nk is not None,
                                               chunk_size=args.chunk_size,
                                               max_iter=args.refractor_max_iter or None,
                                               identity_max_iter=args.identity_max_iter or None,
                                               tolerance=args.refractor_tolerance)
            if not all(np.isfinite(value).all() for value in solution.values()):
                raise RuntimeError(f"Nonfinite refractor solution: NK={nk}, epsilon={epsilon}")
            # Diagnose feasibility without changing the original stopping rule.
            source_mass = np.empty(nk)
            log_target_mass = np.full(target_nk, -np.inf)
            for start in range(0, nk, 512):
                stop = min(start + 512, nk)
                log_plan = (np.log(p[start:stop, None]) + np.log(q[None, :])
                            + (solution["f_raw"][start:stop, None]
                               + solution["g_raw"][None, :]
                               - cost_matrix_chunk(x[start:stop], y))/epsilon)
                source_mass[start:stop] = np.exp(logsumexp(log_plan, axis=1))
                log_target_mass = np.logaddexp(log_target_mass, logsumexp(log_plan, axis=0))
            diagnostic = {"nk": nk, "target_nk": target_nk, "epsilon": epsilon,
                          "source_marginal_l1": float(np.abs(source_mass-p).sum()),
                          "target_marginal_l1": float(np.abs(np.exp(log_target_mass)-q).sum())}
            diagnostic.update({name: solution[name] for name in (
                "final_iterations", "final_change", "source_identity_iterations",
                "source_identity_change", "target_identity_iterations", "target_identity_change")})
            (args.output / f"residual_NK{nk}{target_tag}_eps{epsilon:.17g}.json").write_text(
                json.dumps(diagnostic, indent=2), encoding="utf-8")
            artifacts = build_refraction_artifacts(
                x, y, p, q, solution["f_raw"], solution["g_raw"],
                solution["f"], solution["g"], P_refraction_patch, Q_refraction_patch,
                grid_res=16,
            )
            cache[key] = (solution, artifacts, y, p, q)
            np.savez(checkpoint, total_cost=solution["total_cost"],
                     y_push_3d=artifacts["y_push_3d"], y=y, p=p, q=q)
        return cache[key]

    def cost_row(experiment, nk, epsilon):
        if (experiment, nk, epsilon) in imported:
            row = imported[experiment, nk, epsilon]
            if not row["pf_epsilon"] or float(row["pf_epsilon"]) != args.cost_epsilon:
                raise ValueError("Imported experiments 1-3 must contain matching PF EOT evaluation epsilon")
            rows.append(row)
            write_costs()
            return float(row["cost"])
        _, artifacts, y, p, q = solve(nk, epsilon)
        target_tag = f"_targetNK{args.target_nk}_self" if args.target_nk is not None else ""
        target_tag += solver_tag
        if args.pf_tolerance != 1e-9:
            target_tag += f"_pftol{args.pf_tolerance:.17g}"
        checkpoint = args.output / f"pf_NK{nk}{target_tag}_eps{epsilon:.17g}_eval{args.cost_epsilon:.17g}.json"
        if checkpoint.exists():
            value, info = json.loads(checkpoint.read_text(encoding="utf-8"))
        else:
            value, info = entropic_ot_cost(artifacts["y_push_3d"], y[q > 0], p[p > 0], q[q > 0],
                                           epsilon=args.cost_epsilon, return_info=True,
                                           chunk_size=args.chunk_size, tolerance=args.pf_tolerance,
                                           max_iter=args.pf_max_iter, max_cache_entries=args.pf_cache_entries)
            if info["converged"]:
                checkpoint.write_text(json.dumps([value, info]), encoding="utf-8")
        if not info["converged"] or not np.isfinite(value):
            raise RuntimeError(f"PF EOT evaluation failed: NK={nk}, epsilon={epsilon}: {info}")
        row = dict(experiment=experiment, nk=nk, refractor_epsilon=epsilon,
                   pf_epsilon=args.cost_epsilon, cost=value,
                   converged=info["converged"], iterations=info["iterations"])
        rows.append(row)
        write_costs()
        return row["cost"]

    def save(fig, name):
        fig.tight_layout()
        fig.savefig(args.output / f"{name}.png", dpi=180)
        plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    specs = [
        (1, args.nk, [joint_epsilon(n) for n in args.nk], "Joint movement", r"$N_k$"),
        (2, args.nk, [args.fixed_epsilon]*len(args.nk),
         f"Fixed epsilon = {args.fixed_epsilon:.6g}", r"$N_k$"),
        (3, [args.fixed_nk]*len(args.epsilon), args.epsilon,
         f"Fixed NK = {args.fixed_nk}", r"Refractor $\epsilon$"),
    ]
    for ax, (experiment, ns, es, title, xlabel) in zip(axes, specs):
        if args.target_nk is not None:
            title += f"\nTarget fixed at {args.target_nk} points"
            xlabel = r"Source $N_k$" if experiment != 3 else xlabel
        if solver_tag:
            title += f"\nTransport cap: {args.refractor_max_iter or 'none'}; self cap: {args.identity_max_iter or 'none'}"
        if experiment not in args.experiments:
            ax.set_visible(False)
            continue
        costs = [cost_row(experiment, n, e) for n, e in zip(ns, es)]
        ax.plot(es if experiment == 3 else ns, costs, "o-")
        ax.set(title=f"{experiment}. {title}", xlabel=xlabel,
               ylabel=rf"$EOT_{{{args.cost_epsilon:.6g}}}(PF,T)$ (squared Euclidean, 3D)")
        ax.grid(alpha=0.3)
        ax.set_xscale("log")
        single, single_ax = plt.subplots(figsize=(7, 5))
        single_ax.plot(es if experiment == 3 else ns, costs, "o-")
        single_ax.set(title=ax.get_title(), xlabel=xlabel,
                      ylabel=ax.get_ylabel(), xscale="log")
        single_ax.grid(alpha=0.3)
        save(single, f"experiment_{experiment}")
    save(fig, "experiments_1_2_3")
    if 4 not in args.experiments:
        write_costs()
        print(f"Saved selected experiments to {args.output}", flush=True)
        return

    # A joint-schedule panel plus the full independent NK x refractor-epsilon sweep.
    schedules = [("Joint epsilon(NK)", None)] + [(f"Refractor epsilon = {e:.6g}", e) for e in args.epsilon]
    if args.joint_only:
        schedules = schedules[:1]
    columns = min(3, len(schedules))
    nrows = int(np.ceil(len(schedules) / columns))
    fig, axes = plt.subplots(nrows, columns, figsize=(10*columns, 6*nrows), squeeze=False)
    ordering = []
    pf_potentials = {}
    pf_nk = [nk for nk in args.nk if nk not in args.skip_pf_nk]
    for ax, (title, fixed_e) in zip(axes.flat, schedules):
        ax.set_prop_cycle(color=plt.cm.viridis(np.linspace(0, 1, len(args.pf_epsilon))))
        for pf_e in sorted(args.pf_epsilon, reverse=args.pf_warm_start):
            costs = []
            for nk in pf_nk:
                e = joint_epsilon(nk) if fixed_e is None else fixed_e
                if (nk, e, pf_e) in imported_pf:
                    row = imported_pf[nk, e, pf_e]
                    rows.append(row)
                    costs.append(float(row["cost"]))
                    continue
                _, artifacts, y, p, q = solve(nk, e)
                target_tag = f"_targetNK{args.target_nk}_self" if args.target_nk is not None else ""
                target_tag += solver_tag
                if args.pf_tolerance != 1e-9:
                    target_tag += f"_pftol{args.pf_tolerance:.17g}"
                checkpoint = args.output / f"pf_NK{nk}{target_tag}_eps{e:.17g}_eval{pf_e:.17g}.json"
                if checkpoint.exists():
                    value, info = json.loads(checkpoint.read_text(encoding="utf-8"))
                else:
                    value, info = entropic_ot_cost(
                        artifacts["y_push_3d"], y[q > 0], p[p > 0], q[q > 0],
                        epsilon=pf_e, return_info=True, max_iter=args.pf_max_iter,
                        chunk_size=args.chunk_size, tolerance=args.pf_tolerance,
                        initial_potentials=pf_potentials.get((nk, e)),
                        return_potentials=args.pf_warm_start,
                        max_cache_entries=args.pf_cache_entries,
                    )
                    if args.pf_warm_start:
                        pf_potentials[nk, e] = info.pop("potentials")
                    if info["converged"]:
                        checkpoint.write_text(json.dumps([value, info]), encoding="utf-8")
                if not info["converged"]:
                    raise RuntimeError(f"PF cost failed to converge: NK={nk}, epsilon={e}, PF={pf_e}: {info}")
                costs.append(value)
                rows.append(dict(experiment=4, nk=nk, refractor_epsilon=e,
                                 pf_epsilon=pf_e, cost=value, **{k: info[k] for k in ("converged", "iterations")}))
                print(f"  PF NK={nk}, eps={pf_e:.6g}: {value:.9g}", flush=True)
            ax.plot(pf_nk, costs, "o-", label=f"PF eps={pf_e:.6g}")
            ordering.append(dict(schedule=title, pf_epsilon=pf_e,
                                 nk_in_ascending_cost_order=[pf_nk[i] for i in np.argsort(costs)],
                                 decreases_with_nk=bool(np.all(np.diff(costs) <= 0))))
        ax.set(title=title, xlabel=r"$N_k$", ylabel=r"$OT_\epsilon(PF,T)$ (squared Euclidean, 3D)")
        ax.set_xscale("log")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7, loc="upper left", bbox_to_anchor=(1.02, 1))
    for ax in list(axes.flat)[len(schedules):]:
        ax.set_visible(False)
    save(fig, "experiment_4_pushforward")
    write_costs()
    (args.output / "ordering.json").write_text(json.dumps(ordering, indent=2), encoding="utf-8")
    (args.output / "parameters.json").write_text(json.dumps(vars(args), default=str, indent=2), encoding="utf-8")
    print(f"Saved plots, costs, parameters and ordering to {args.output}", flush=True)


if __name__ == "__main__":
    main()
