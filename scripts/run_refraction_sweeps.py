"""Four parameter sweeps using the original uniform refractor pipeline.

Run from the repository root: python scripts/run_refraction_sweeps.py
"""
try:
    from result_paths import artifact_path, artifact_glob, resolve_path, dated_run_dir, ensure_readme, self_cost_dir
except ModuleNotFoundError:
    from scripts.result_paths import artifact_path, artifact_glob, resolve_path, dated_run_dir, ensure_readme, self_cost_dir
import argparse
import sys
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
from refracter.optics import entropic_surface_normals, refract_using_normals
import refracter.sinkhorn as sinkhorn_module
from atomic_io import atomic_replace, atomic_write_json
from shared_optical import problem_key, get_or_solve
from refracter.cost_cache import bounded_cost_cache


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
    parser.add_argument("--kappa", type=float, default=0.6)
    for name, bounds in (("source-theta", SRC_THETA), ("source-phi", SRC_PHI),
                         ("target-theta", TGT_THETA), ("target-phi", TGT_PHI)):
        parser.add_argument("--" + name, nargs=2, type=float,
                            default=np.rad2deg(bounds).tolist(), metavar=("MIN", "MAX"),
                            help="Patch bounds in degrees; theta is azimuth, phi is polar")
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
    parser.add_argument("--optical-marginal-tolerance", type=float, default=1e-7,
                        help="Required L1 feasibility of optical transport and both self solves")
    parser.add_argument("--pf-tolerance", type=float, default=1e-9)
    parser.add_argument("--chunk-size", type=int, default=512)
    parser.add_argument("--pf-cache-entries", type=int, default=4_000_000)
    parser.add_argument("--workers", type=int, default=1, help="Threads for independent SciPy reductions")
    parser.add_argument('--shared-optical-cache',type=Path,
                        help='Reuse identical optical potentials between hard and entropic consumers')
    parser.add_argument('--cost-cache-mib',type=int,default=0,
                        help='Maximum retained optical cost blocks per solve (execution setting)')
    parser.add_argument("--pushforward", choices=["hard", "entropic"], default="hard",
                        help="Hard c-transform or entropic-surface normals with Snell refraction")
    parser.add_argument("--pf-warm-start", action="store_true",
                        help="Reuse PF potentials from the preceding larger epsilon")
    parser.add_argument("--output", type=Path, default=Path(REPO_ROOT) / "results/runs/standard_fixed_target_oct_22_26")
    args = parser.parse_args()
    if args.cost_cache_mib<0:
        parser.error('cost-cache-mib must be nonnegative')
    args.solver_version = 'untruncated_marginal_v2'
    if not np.isfinite(args.kappa) or not 0 < args.kappa < 1:
        parser.error("Refractor kappa must be finite and in (0, 1)")
    for name in ("source_theta", "source_phi", "target_theta", "target_phi"):
        lo, hi = getattr(args, name)
        limit = 180 if name.endswith("phi") else 360
        if not np.isfinite([lo, hi]).all() or not 0 <= lo < hi <= limit:
            parser.error(f"{name} must satisfy 0 <= MIN < MAX <= {limit}")
    def patch_density(theta_bounds, phi_bounds):
        def density(points):
            theta = np.mod(np.rad2deg(np.arctan2(points[:, 1], points[:, 0])), 360)
            phi = np.rad2deg(np.arccos(np.clip(points[:, 2], -1, 1)))
            theta_inside = ((theta >= theta_bounds[0]-1e-10) & (theta <= theta_bounds[1]+1e-10))
            theta_inside |= np.hypot(points[:,0],points[:,1]) <= 1e-14
            return (theta_inside
                    & (phi >= phi_bounds[0]-1e-10) & (phi <= phi_bounds[1]+1e-10)).astype(float)
        return density
    source_density = patch_density(args.source_theta, args.source_phi)
    target_density = patch_density(args.target_theta, args.target_phi)
    args.output = resolve_path(args.output)
    if not any(value == "--output" or value.startswith("--output=") for value in sys.argv[1:]):
        label = "uniform_refraction_hard_map" if args.pushforward == "hard" else "smooth_refraction"
        if args.extended:
            label += "_extended"
        elif args.expanded:
            label += "_expanded"
        args.output = dated_run_dir("runs", label)
    if args.import_costs:
        args.import_costs = resolve_path(args.import_costs)
    if min(args.refractor_max_iter, args.identity_max_iter) < 0:
        parser.error("Iteration caps must be nonnegative (0 means uncapped)")
    if args.chunk_size < 1 or not all(np.isfinite(v) and v > 0 for v in (args.refractor_tolerance, args.pf_tolerance, args.optical_marginal_tolerance)):
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
    parameter_path = artifact_path(args.output, "parameters.json")
    if parameter_path.exists():
        previous = json.loads(parameter_path.read_text(encoding="utf-8"))
        defaults = dict(kappa=0.6, source_theta=np.rad2deg(SRC_THETA).tolist(),
                        source_phi=np.rad2deg(SRC_PHI).tolist(),
                        target_theta=np.rad2deg(TGT_THETA).tolist(),
                        target_phi=np.rad2deg(TGT_PHI).tolist())
        for name, default in defaults.items():
            if not np.allclose(previous.get(name, default), getattr(args, name), rtol=0, atol=1e-12):
                parser.error(f"Existing output has different {name}; use a separate output folder")
        for name in ('solver_version','optical_marginal_tolerance','refractor_max_iter',
                     'identity_max_iter','refractor_tolerance','pf_tolerance','pf_max_iter',
                     'target_nk','fixed_nk','fixed_epsilon','cost_epsilon','pushforward',
                     'nk','epsilon','pf_epsilon','experiments'):
            if previous.get(name) != getattr(args,name):
                parser.error(f'Existing output has different {name}; use a separate output folder')
    args.self_cost_mode = "separate_marginals_v1" if args.target_nk is not None else "legacy"
    args.cost_metric = "pushforward_target_eot_squared_euclidean"
    atomic_write_json(artifact_path(args.output, "parameters.json", for_write=True), vars(args))
    set_kappa(args.kappa)
    cache = {}
    pf_potentials = {}
    solver_tag = "_untruncated_v2"
    if (args.refractor_max_iter, args.identity_max_iter, args.refractor_tolerance) != (17, 17, 1e-5):
        solver_tag += f"_iter{args.refractor_max_iter}_id{args.identity_max_iter}_tol{args.refractor_tolerance:.17g}"
    if args.pushforward == "entropic":
        solver_tag += "_entropic_v1"
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

    def checkpoint_tag(nk, *, evaluation=False):
        tag = f"_targetNK{args.target_nk or nk}_self" if args.target_nk is not None else ""
        tag += solver_tag
        if evaluation and args.pf_tolerance != 1e-9:
            tag += f"_pftol{args.pf_tolerance:.17g}"
        return tag

    def write_costs():
        temporary = artifact_path(args.output, "costs.tmp.csv", for_write=True)
        with temporary.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)
        atomic_replace(temporary, artifact_path(args.output, "costs.csv", for_write=True))

    def solve(nk, epsilon):
        key = (nk, epsilon)
        target_nk = args.target_nk or nk
        target_tag = checkpoint_tag(nk)
        if key not in cache:
            checkpoint = artifact_path(args.output, f"solve_NK{nk}{target_tag}_eps{epsilon:.17g}.npz")
            if checkpoint.exists():
                residual = artifact_path(args.output, f"residual_NK{nk}{target_tag}_eps{epsilon:.17g}.json")
                diagnostic = json.loads(residual.read_text(encoding='utf-8'))
                errors = np.array([diagnostic[name] for name in ('source_marginal_l1','target_marginal_l1',
                    'source_identity_marginal_l1','target_identity_marginal_l1')])
                if not np.isfinite(errors).all() or errors.max() > args.optical_marginal_tolerance:
                    raise RuntimeError('Saved optical checkpoint fails marginal feasibility')
                with np.load(checkpoint) as data:
                    cache[key] = ({"total_cost": float(data["total_cost"])},
                                  {"y_push_3d": data["y_push_3d"]},
                                  data["y"], data["p"], data["q"])
                return cache[key]
            print(f"Refractor NK={nk}, epsilon={epsilon:.8g}", flush=True)
            x = gen_spherical_patch(nk, *np.deg2rad(args.source_theta), *np.deg2rad(args.source_phi), skip=0)
            y = gen_spherical_patch(target_nk, *np.deg2rad(args.target_theta), *np.deg2rad(args.target_phi), skip=0)
            p, q = source_density(x), target_density(y)
            p, q = p/p.sum(), q/q.sum()
            def calculate_optics(xx,yy):
                return run_sinkhorn_divergence(xx, yy, p, q, verbose=False, epsilon=epsilon,
                                               separate_marginals=args.target_nk is not None,
                                               chunk_size=args.chunk_size,
                                               max_iter=args.refractor_max_iter or None,
                                               identity_max_iter=args.identity_max_iter or None,
                                               tolerance=args.refractor_tolerance,
                                               marginal_tolerance=args.optical_marginal_tolerance)
            def optical_solve():
                if not args.cost_cache_mib:
                    return calculate_optics(x,y)
                with bounded_cost_cache(x,y,args.cost_cache_mib*1024**2,args.chunk_size) as costs:
                    result=calculate_optics(costs.x,costs.y)
                    execution=artifact_path(args.output,f'cost_cache_NK{nk}_eps{epsilon:.17g}.json',for_write=True)
                    atomic_write_json(execution,costs.stats())
                    print(f'Cost cache NK={nk}: {costs.stats()}',flush=True)
                    return result
            if args.shared_optical_cache:
                settings=dict(solver_version=args.solver_version,kappa=args.kappa,epsilon=epsilon,
                    separate_marginals=args.target_nk is not None,chunk_size=args.chunk_size,
                    max_iter=args.refractor_max_iter,identity_max_iter=args.identity_max_iter,
                    tolerance=args.refractor_tolerance,marginal_tolerance=args.optical_marginal_tolerance)
                identity=problem_key((x,y,p,q),settings)
                solution=get_or_solve(args.shared_optical_cache,identity,optical_solve)
            else:
                solution=optical_solve()
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
            diagnostic.update({name:solution[name] for name in (
                'source_identity_marginal_l1','target_identity_marginal_l1')})
            if max(diagnostic['source_marginal_l1'],diagnostic['target_marginal_l1']) > args.optical_marginal_tolerance:
                raise RuntimeError(f'Optical marginal feasibility failed: {diagnostic}')
            artifacts = build_refraction_artifacts(
                x, y, p, q, solution["f_raw"], solution["g_raw"],
                solution["f"], solution["g"], source_density, target_density,
                grid_res=16, kappa=args.kappa,
            )
            entropic_arrays = {}
            if args.pushforward == "entropic":
                normals, log_radius = entropic_surface_normals(
                    x, y, q, solution["g_raw"], epsilon, args.kappa, args.chunk_size)
                transmitted = refract_using_normals(x[p > 0], normals[p > 0], args.kappa)
                if not np.isfinite(transmitted).all() or not np.allclose(
                        np.linalg.norm(transmitted, axis=1), 1, atol=1e-10, rtol=0):
                    raise RuntimeError("Nonfinite or nonunit entropic transmitted directions")
                entropic_arrays = dict(x=x, f_raw=solution["f_raw"], g_raw=solution["g_raw"],
                                       surface_normals=normals, log_radius=log_radius,
                                       hard_y_push_3d=artifacts["y_push_3d"].copy())
                artifacts["y_push_3d"] = transmitted
                diagnostic.update(pushforward="entropic", normal_dot_error=float(
                    np.max(np.abs(np.sum(x*normals, axis=1)-1))),
                    transmitted_norm_error=float(np.max(np.abs(np.linalg.norm(transmitted, axis=1)-1))))
            (artifact_path(args.output, f"residual_NK{nk}{target_tag}_eps{epsilon:.17g}.json", for_write=True)).write_text(
                json.dumps(diagnostic, indent=2), encoding="utf-8")
            cache[key] = (solution, artifacts, y, p, q)
            np.savez(artifact_path(args.output, checkpoint.name, for_write=True), total_cost=solution["total_cost"],
                     y_push_3d=artifacts["y_push_3d"], y=y, p=p, q=q, **entropic_arrays)
        return cache[key]

    def evaluate_cost(nk, epsilon, pf_epsilon, *, warm_start=False):
        """One PF EOT evaluation/cache path shared by all four experiments."""
        _, artifacts, y, p, q = solve(nk, epsilon)
        tag = checkpoint_tag(nk, evaluation=True)
        checkpoint = artifact_path(args.output, f"pf_NK{nk}{tag}_eps{epsilon:.17g}_eval{pf_epsilon:.17g}.json")
        if checkpoint.exists():
            value, info = json.loads(checkpoint.read_text(encoding="utf-8"))
        else:
            value, info = entropic_ot_cost(
                artifacts["y_push_3d"], y[q > 0], p[p > 0], q[q > 0],
                epsilon=pf_epsilon, return_info=True, chunk_size=args.chunk_size,
                tolerance=args.pf_tolerance, max_iter=args.pf_max_iter,
                max_cache_entries=args.pf_cache_entries,
                initial_potentials=pf_potentials.get((nk, epsilon)) if warm_start else None,
                return_potentials=warm_start)
            if warm_start:
                pf_potentials[nk, epsilon] = info.pop("potentials")
            if info["converged"] and np.isfinite(value):
                artifact_path(args.output, checkpoint.name, for_write=True).write_text(json.dumps([value, info]), encoding="utf-8")
        if not info["converged"] or not np.isfinite(value):
            raise RuntimeError(f"PF EOT evaluation failed: NK={nk}, refractor={epsilon}, PF={pf_epsilon}: {info}")
        return value, info

    def cost_row(experiment, nk, epsilon):
        if (experiment, nk, epsilon) in imported:
            row = imported[experiment, nk, epsilon]
            if not row["pf_epsilon"] or float(row["pf_epsilon"]) != args.cost_epsilon:
                raise ValueError("Imported experiments 1-3 must contain matching PF EOT evaluation epsilon")
            rows.append(row)
            write_costs()
            return float(row["cost"])
        value, info = evaluate_cost(nk, epsilon, args.cost_epsilon)
        row = dict(experiment=experiment, nk=nk, refractor_epsilon=epsilon,
                   pf_epsilon=args.cost_epsilon, cost=value,
                   converged=info["converged"], iterations=info["iterations"])
        rows.append(row)
        write_costs()
        return row["cost"]

    def save(fig, name):
        fig.tight_layout()
        fig.savefig(artifact_path(args.output, f"{name}.png", for_write=True), dpi=180)
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
        title += f"\nPush-forward: {args.pushforward}"
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
                value, info = evaluate_cost(nk, e, pf_e, warm_start=args.pf_warm_start)
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
    (artifact_path(args.output, "ordering.json", for_write=True)).write_text(json.dumps(ordering, indent=2), encoding="utf-8")
    atomic_write_json(artifact_path(args.output, "parameters.json", for_write=True), vars(args))
    print(f"Saved plots, costs, parameters and ordering to {args.output}", flush=True)


if __name__ == "__main__":
    main()
