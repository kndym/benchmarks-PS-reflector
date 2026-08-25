# Shared implementation for the point-source far-field refraction benchmark.
# The command-line scripts choose inputs and output names; the numerical work
# stays here so the single case and density sweep use the same pipeline.

from __future__ import annotations

import os
import time
from typing import Callable, Dict, Optional

import numpy as np

from .build import c_transform_fc, c_transform_gc
from .cost import cost_matrix_chunk, set_kappa
from .distributions import (
    P_refraction_patch,
    Q_refraction_patch,
    make_patch_cross,
    make_patch_donut,
    make_patch_gaussian,
    make_patch_uniform,
    stereo_north,
)
from .sinkhorn import (
    sinkhorn_identity_f_step,
    sinkhorn_identity_g_step,
    sinkhorn_step,
)


# The cost is c(x,y) = -log(1 - κ·(x·y)); the paper's refractor
# experiments use kappa = 0.6.
KAPPA = 0.6
# Angular bounds for the source and target spherical patches.
SRC_THETA = (np.pi / 12, np.pi / 3)
SRC_PHI = (np.pi / 12, np.pi / 4)
TGT_THETA = (np.pi / 3, 5 * np.pi / 12)
TGT_PHI = (np.pi / 10, np.pi / 5)

DENSITY_FACTORIES = {
    "uniform": make_patch_uniform,
    "gaussian": make_patch_gaussian,
    "donut": make_patch_donut,
    "cross": make_patch_cross,
}
# The all-pairs script runs every ordered pair of these four density shapes.
DENSITY_NAMES = ("uniform", "gaussian", "donut", "cross")


# Return one low-discrepancy coordinate for the spherical point cloud.
def _halton(index: int, base: int) -> float:
    result, fraction, remainder = 0.0, 1.0, index
    while remainder > 0:
        fraction /= base
        result += fraction * (remainder % base)
        remainder //= base
    return result


def gen_spherical_patch(
    n: int,
    theta_min: float,
    theta_max: float,
    phi_min: float,
    phi_max: float,
    *,
    base2: int = 2,
    base3: int = 3,
    skip: int = 0,
) -> np.ndarray:
    # Map Halton points to spherical angles, then convert them to unit vectors.
    points = []
    index = skip
    while len(points) < n:
        u1 = _halton(index, base2)
        u2 = _halton(index, base3)
        phi = phi_min + u2 * (phi_max - phi_min)
        theta = theta_min + u1 * (theta_max - theta_min)
        points.append(
            [
                np.sin(phi) * np.cos(theta),
                np.sin(phi) * np.sin(theta),
                np.cos(phi),
            ]
        )
        index += 1
    return np.asarray(points, dtype=np.float64)


def make_patch_density(name: str, *, source: bool) -> Callable[[np.ndarray], np.ndarray]:
    # Select a density and apply the source or target patch geometry.
    if name not in DENSITY_FACTORIES:
        choices = ", ".join(sorted(DENSITY_FACTORIES))
        raise ValueError(f"Unknown density {name!r}; choose from: {choices}")
    theta = SRC_THETA if source else TGT_THETA
    phi = SRC_PHI if source else TGT_PHI
    return DENSITY_FACTORIES[name](*theta, *phi)


def _run_identity_f(
    x: np.ndarray,
    y: np.ndarray,
    logp: np.ndarray,
    p: np.ndarray,
    k_final: int,
    id_step: int,
    chunk_size: int,
    cap_iter: int,
    verbose: bool,
) -> np.ndarray:
    # Solve the source self-transport problem used in the Sinkhorn divergence.
    f_id = np.zeros(len(x), dtype=np.float64)
    regvar = 1
    while regvar < k_final:
        f_id_new, _ = sinkhorn_identity_f_step(x, y, logp, f_id, regvar, chunk_size)
        f_id = np.where(p > 0, f_id_new, f_id)
        regvar += id_step

    maxdif = 1.0 + 1e-5
    iterations = 0
    while maxdif > 1e-5 and iterations < cap_iter:
        f_id_new, maxdif = sinkhorn_identity_f_step(
            x, y, logp, f_id, k_final, chunk_size
        )
        f_id = np.where(p > 0, f_id_new, f_id)
        iterations += 1
    if verbose:
        print(f"  Identity F: {iterations} final iterations, maxdif={maxdif:.3e}")
    return f_id


def _run_identity_g(
    x: np.ndarray,
    y: np.ndarray,
    logq: np.ndarray,
    q: np.ndarray,
    k_final: int,
    id_step: int,
    chunk_size: int,
    cap_iter: int,
    verbose: bool,
) -> np.ndarray:
    # Solve the target self-transport problem used in the Sinkhorn divergence.
    g_id = np.zeros(len(y), dtype=np.float64)
    regvar = 1
    while regvar < k_final:
        g_id_new, _ = sinkhorn_identity_g_step(x, y, logq, g_id, regvar, chunk_size)
        g_id = np.where(q > 0, g_id_new, g_id)
        regvar += id_step

    maxdif = 1.0 + 1e-5
    iterations = 0
    while maxdif > 1e-5 and iterations < cap_iter:
        g_id_new, maxdif = sinkhorn_identity_g_step(
            x, y, logq, g_id, k_final, chunk_size
        )
        g_id = np.where(q > 0, g_id_new, g_id)
        iterations += 1
    if verbose:
        print(f"  Identity G: {iterations} final iterations, maxdif={maxdif:.3e}")
    return g_id


def run_refraction_case(
    source_density: Callable[[np.ndarray], np.ndarray],
    target_density: Callable[[np.ndarray], np.ndarray],
    *,
    nk: int = 1600,
    chunk_size: int = 512,
    use_identity: bool = False,
    verbose: bool = True,
) -> Dict[str, np.ndarray]:
    # Run one source-to-target case and return arrays used by the notebook.
    set_kappa(KAPPA)
    t_start = time.time()

    x = gen_spherical_patch(nk, *SRC_THETA, *SRC_PHI)
    y = gen_spherical_patch(nk, *TGT_THETA, *TGT_PHI)
    p_raw = np.asarray(source_density(x), dtype=np.float64)
    q_raw = np.asarray(target_density(y), dtype=np.float64)
    if p_raw.sum() <= 0 or q_raw.sum() <= 0:
        raise ValueError("Source and target densities must both have positive mass")

    p = p_raw / p_raw.sum()
    q = q_raw / q_raw.sum()
    # Log-domain Sinkhorn needs -inf for points outside a density's support.
    logp = np.full_like(p, -np.inf)
    logq = np.full_like(q, -np.inf)
    np.log(p, out=logp, where=p > 0)
    np.log(q, out=logq, where=q > 0)

    # Match the reference regularization schedule: k_final = 8 floor(sqrt(NK)).
    multiplier = 8
    k_final = multiplier * int(np.floor(np.sqrt(nk)))
    id_step = int(np.floor(np.sqrt(k_final)))
    cap_iter = 16
    cap_thr = 1e-5

    if verbose:
        print(f"Refraction case: NK={nk}, k_final={k_final}")
        print(f"  Source support: {(p > 0).sum()} / {nk}")
        print(f"  Target support: {(q > 0).sum()} / {nk}")

    # Main entropic OT solve between the source and target point clouds.
    f = np.zeros(nk, dtype=np.float64)
    g = np.zeros(nk, dtype=np.float64)
    maxdif = cap_thr + 1.0
    iterations = 0
    while maxdif > cap_thr and iterations < cap_iter:
        f, g, maxdif = sinkhorn_step(x, y, logp, logq, f, g, k_final, chunk_size)
        iterations += 1
    if verbose:
        print(f"  Sinkhorn: {iterations} iterations, maxdif={maxdif:.3e}")

    # Keep the raw OT potentials for the c-transforms and push-forward map.
    f_raw = f.copy()
    g_raw = g.copy()
    f_id = np.zeros(nk, dtype=np.float64)
    g_id = np.zeros(nk, dtype=np.float64)

    if use_identity:
        # Subtract source and target self-costs for the full Sinkhorn divergence.
        f_id = _run_identity_f(
            x, y, logp, p, k_final, id_step, chunk_size, cap_iter, verbose
        )
        g_id = _run_identity_g(
            x, y, logq, q, k_final, id_step, chunk_size, cap_iter, verbose
        )

        f = np.where(p > 0, f, 0.0)
        g = np.where(q > 0, g, 0.0)
        max_f_id = float(np.max(f_id))
        f_id -= max_f_id
        g_id += max_f_id
        max_f = float(np.max(f))
        f -= max_f
        g += max_f
        f -= f_id
        g -= g_id

    # c-transforms verify the dual potentials and define the surface geometry.
    mask_p = p > 0
    mask_q = q > 0
    gc = c_transform_gc(
        x, y, np.where(mask_p, f_raw, -np.inf), chunk_size=chunk_size
    )
    fc = c_transform_fc(
        x, y, np.where(mask_q, g_raw, -np.inf), chunk_size=chunk_size
    )
    # The refractor radius is R = exp(f); the surface point is 2 R x.
    R = np.exp(f)
    Ref = 2.0 * x * R[:, None]
    Refc = 2.0 * x * np.exp(fc)[:, None]

    # Evaluate both densities on a regular stereographic grid for plotting.
    grid_res = 256
    grid_side = np.linspace(-0.6, 0.6, grid_res)
    uu, vv = np.meshgrid(grid_side, grid_side, indexing="ij")
    n2 = uu**2 + vv**2
    denom = 1.0 + n2
    grid = np.stack(
        [2 * uu / denom, 2 * vv / denom, (1 - n2) / denom], axis=-1
    ).reshape(-1, 3)
    X_MeshGrid = np.asarray(source_density(grid)).reshape(grid_res, grid_res)
    Y_MeshGrid = np.asarray(target_density(grid)).reshape(grid_res, grid_res)

    u_y, v_y = stereo_north(y[mask_q])
    Y_projected = np.column_stack([u_y, v_y, q[mask_q]])

    # Map each supported source point to the target minimizing C(x,y) - g(y).
    src_idx = np.where(mask_p)[0]
    tgt_idx = np.where(mask_q)[0]
    pushed_y = []
    for start in range(0, len(src_idx), chunk_size):
        rows = src_idx[start : start + chunk_size]
        costs = cost_matrix_chunk(x[rows], y[tgt_idx])
        j_star = np.argmin(costs - g_raw[tgt_idx][None, :], axis=1)
        pushed_y.append(y[tgt_idx[j_star]])
    pushed_y = np.vstack(pushed_y)
    u_push, v_push = stereo_north(pushed_y)
    Y_Pushed_projected = np.column_stack(
        [u_push, v_push, np.asarray(target_density(pushed_y))]
    )

    # Discrete OT objective, using corrected potentials when identity terms exist.
    total_cost = float(
        np.sum(p[mask_p] * f[mask_p]) + np.sum(q[mask_q] * g[mask_q])
    )
    if verbose:
        print(f"  Total cost: {total_cost:.6e}")
        print(f"  Elapsed: {time.time() - t_start:.1f}s")

    return {
        "x": x,
        "y": y,
        "p": p,
        "q": q,
        "f_raw": f_raw,
        "g_raw": g_raw,
        "f_id": f_id,
        "g_id": g_id,
        "f": f,
        "g": g,
        "R": R,
        "Ref": Ref,
        "gc": gc,
        "fc": fc,
        "Refc": Refc,
        "X_MeshGrid": X_MeshGrid,
        "Y_MeshGrid": Y_MeshGrid,
        "grid_side": grid_side,
        "Y_projected": Y_projected,
        "Y_Pushed_projected": Y_Pushed_projected,
        "y_push_3d": pushed_y,
        "kappa": np.float64(KAPPA),
        "total_cost": np.float64(total_cost),
    }


def save_refraction_result(
    result: Dict[str, np.ndarray],
    path: str,
    *,
    source_density: Optional[str] = None,
    target_density: Optional[str] = None,
) -> None:
    # Save the stable NPZ schema consumed by Benchmark_refraction.ipynb.
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    payload = {
        key: result[key]
        for key in (
            "x", "y", "p", "q", "f_raw", "g_raw", "f_id", "g_id",
            "f", "g", "R", "Ref", "gc", "fc", "Refc", "X_MeshGrid",
            "Y_MeshGrid", "grid_side", "Y_projected", "Y_Pushed_projected",
            "y_push_3d", "kappa",
        )
    }
    payload["x_s"] = np.zeros((0, 3), dtype=np.float64)
    payload["y_s"] = np.zeros((0, 3), dtype=np.float64)
    if source_density is not None:
        payload["src_density"] = source_density
    if target_density is not None:
        payload["tgt_density"] = target_density
    np.savez(path, **payload)


def save_surface_figure(result: Dict[str, np.ndarray], path: str) -> None:
    # Save a quick 3-D view of the computed refractor surface.
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection="3d")
    scatter = ax.scatter(
        result["Ref"][:, 0],
        result["Ref"][:, 1],
        result["Ref"][:, 2],
        c=result["R"],
        cmap="viridis",
        s=15,
        vmin=result["R"].min(),
        vmax=result["R"].max(),
        depthshade=True,
    )
    plt.colorbar(scatter, ax=ax, label="R (refractor radius)", shrink=0.65)
    ax.set_title(f"Refractor Surface - Python (NK={len(result['x'])}, kappa=0.6)")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
