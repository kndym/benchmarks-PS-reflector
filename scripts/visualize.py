#!/usr/bin/env python3
"""Reflector benchmark — static matplotlib visualization.

Saves PNG figures to the output directory.

Usage:
    python visualize.py --output-dir <path> [--benchmark <name>] [--nk <int>]

Called automatically by Benchmark.ipynb (Step 6) but can be re-run standalone
at any time to regenerate the figures without re-running the full benchmark.
"""

import os
import sys
import argparse

import numpy as np
import matplotlib

matplotlib.use("Agg")  # non-interactive backend — safe in subprocess and CLI
import matplotlib.pyplot as plt


# ── Argument parsing ──────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description="Visualize reflector benchmark output")
    p.add_argument("--output-dir", required=True,
                   help="Path to the benchmark output directory (contains .txt files)")
    p.add_argument("--benchmark", default="",
                   help="Benchmark header filename (for plot title)")
    p.add_argument("--nk", type=int, default=0,
                   help="Point-cloud size NK (for plot title)")
    return p.parse_args()


# ── Data loaders ──────────────────────────────────────────────────────────────

def load_meshgrid(path):
    """Load a grid file: every non-empty line is a row of space-separated floats."""
    if not os.path.exists(path):
        return None
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.strip():
                rows.append([float(v) for v in line.split()])
    return np.array(rows) if rows else None


def load_vector(path):
    """Load a file whose first line is a header; remaining lines are float rows."""
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        fh.readline()               # skip header line
        rows = []
        for line in fh:
            if line.strip():
                vals = [float(v) for v in line.split()]
                if vals:
                    rows.append(vals)
    return np.array(rows) if rows else None


def load_points(path):
    """Load a file of rows with >= 2 floats (no header)."""
    if not os.path.exists(path):
        return None
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.strip():
                vals = [float(v) for v in line.split()]
                if len(vals) >= 2:
                    rows.append(vals)
    return np.array(rows) if rows else None


# ── Helpers ───────────────────────────────────────────────────────────────────

def no_data(ax, msg):
    ax.text(0.5, 0.5, msg, ha="center", va="center",
            transform=ax.transAxes, color="grey", fontsize=10)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    args      = parse_args()
    out_dir   = args.output_dir
    BENCHMARK = args.benchmark
    NK        = args.nk

    j = lambda name: os.path.join(out_dir, name)

    # Load all output files
    X_mesh = load_meshgrid(j("X_MeshGrid.txt"))
    Y_mesh = load_meshgrid(j("Y_MeshGrid.txt"))
    Y_proj = load_points(j("Y_projected.txt"))
    _yp    = load_points(j("Y_Pushed_projected.txt"))
    Y_push = _yp if _yp is not None else load_points(j("Y_Pushed_projected_owndisc.txt"))
    Ref    = load_vector(j("Ref_MY.txt"))
    R_data = load_vector(j("R_MY.txt"))

    # ── 6-panel main figure ───────────────────────────────────────────────────
    fig = plt.figure(figsize=(20, 13))
    fig.suptitle(
        f"Reflector Benchmark — {BENCHMARK} / NK={NK}",
        fontsize=13, fontweight="bold",
    )

    # Panel 1 — source density
    ax = plt.subplot(2, 3, 1)
    if X_mesh is not None:
        im = ax.imshow(X_mesh, extent=[-0.6, 0.6, -0.6, 0.6],
                       origin="lower", cmap="viridis")
        plt.colorbar(im, ax=ax, label="Density")
        ax.set_title("Source Density", fontweight="bold")
    else:
        no_data(ax, "X_MeshGrid.txt not found")
    ax.set_xlabel("X"); ax.set_ylabel("Y")

    # Panel 2 — destination density
    ax = plt.subplot(2, 3, 2)
    if Y_mesh is not None:
        im = ax.imshow(Y_mesh, extent=[-0.6, 0.6, -0.6, 0.6],
                       origin="lower", cmap="plasma")
        plt.colorbar(im, ax=ax, label="Density")
        ax.set_title("Destination Density", fontweight="bold")
    else:
        no_data(ax, "Y_MeshGrid.txt not found")
    ax.set_xlabel("X"); ax.set_ylabel("Y")

    # Panel 3 — target distribution (original)
    ax = plt.subplot(2, 3, 3)
    if Y_proj is not None and len(Y_proj):
        c  = Y_proj[:, 2] if Y_proj.shape[1] >= 3 else "blue"
        sc = ax.scatter(Y_proj[:, 0], Y_proj[:, 1],
                        c=c, cmap="coolwarm", s=1, alpha=0.6)
        if Y_proj.shape[1] >= 3:
            plt.colorbar(sc, ax=ax, label="Z")
        ax.set_title("Target (original)", fontweight="bold")
        ax.set_aspect("equal"); ax.grid(True, alpha=0.3)
    else:
        no_data(ax, "Y_projected.txt not found")
    ax.set_xlabel("X"); ax.set_ylabel("Y")

    # Panel 4 — pushed / reflected distribution
    ax = plt.subplot(2, 3, 4)
    if Y_push is not None and len(Y_push):
        c  = Y_push[:, 2] if Y_push.shape[1] >= 3 else "red"
        sc = ax.scatter(Y_push[:, 0], Y_push[:, 1],
                        c=c, cmap="coolwarm", s=1, alpha=0.6)
        if Y_push.shape[1] >= 3:
            plt.colorbar(sc, ax=ax, label="Z")
        ax.set_title("Pushed (reflected)", fontweight="bold")
        ax.set_aspect("equal"); ax.grid(True, alpha=0.3)
    else:
        no_data(ax, "Y_Pushed_projected.txt not found")
    ax.set_xlabel("X"); ax.set_ylabel("Y")

    # Panel 5 — reflector surface (3D if available, 2D fallback)
    try:
        from mpl_toolkits.mplot3d import Axes3D
        _has_3d = True
    except Exception:
        _has_3d = False

    if _has_3d:
        ax = plt.subplot(2, 3, 5, projection="3d")
        if Ref is not None and len(Ref) and Ref.shape[1] >= 3:
            idx = np.random.choice(len(Ref), min(len(Ref), 10_000), replace=False)
            ax.scatter(Ref[idx, 0], Ref[idx, 1], Ref[idx, 2],
                       s=0.5, alpha=0.6, c=Ref[idx, 2], cmap="viridis")
            ax.set_title("Reflector Surface (3D)", fontweight="bold")
            ax.set_xlabel("X"); ax.set_ylabel("Y"); ax.set_zlabel("Z")
        else:
            ax.text2D(0.5, 0.5, "Ref_MY.txt not found",
                      ha="center", va="center", transform=ax.transAxes)
    else:
        ax = plt.subplot(2, 3, 5)
        if Ref is not None and len(Ref) and Ref.shape[1] >= 3:
            idx = np.random.choice(len(Ref), min(len(Ref), 10_000), replace=False)
            sc  = ax.scatter(Ref[idx, 0], Ref[idx, 1],
                             c=Ref[idx, 2], cmap="viridis", s=0.5, alpha=0.6)
            plt.colorbar(sc, ax=ax, label="Z")
            ax.set_title("Reflector Surface (2D, Z=colour)", fontweight="bold")
            ax.set_aspect("equal"); ax.grid(True, alpha=0.3)
        else:
            no_data(ax, "Ref_MY.txt not found")
    ax.set_xlabel("X"); ax.set_ylabel("Y")

    # Panel 6 — target vs pushed overlay
    ax = plt.subplot(2, 3, 6)
    if Y_proj is not None and Y_push is not None:
        ax.scatter(Y_proj[:, 0], Y_proj[:, 1],
                   s=1, alpha=0.3, label="Target", color="steelblue")
        ax.scatter(Y_push[:, 0], Y_push[:, 1],
                   s=1, alpha=0.3, label="Pushed", color="tomato")
        ax.legend(markerscale=5)
        ax.set_aspect("equal"); ax.grid(True, alpha=0.3)
        ax.set_title("Target vs Pushed", fontweight="bold")
    else:
        no_data(ax, "comparison data unavailable")
    ax.set_xlabel("X"); ax.set_ylabel("Y")

    plt.tight_layout()
    plt.savefig(j("visualization.png"), dpi=150, bbox_inches="tight")
    print("SAVED:", j("visualization.png"))

    # ── Radius analysis ───────────────────────────────────────────────────────
    if R_data is not None and Ref is not None:
        R_flat = R_data.flatten()
        # Guard against length mismatch between Ref and R_data
        n_common = min(len(Ref), len(R_flat))

        fig, axes = plt.subplots(1, 2, figsize=(16, 6))
        fig.suptitle("Reflector Radius Analysis", fontweight="bold")

        axes[0].hist(R_flat, bins=60, edgecolor="black", alpha=0.7, color="skyblue")
        mean_r = np.mean(R_flat)
        axes[0].axvline(mean_r, color="red", ls="--", lw=2,
                        label=f"Mean = {mean_r:.4f}")
        axes[0].legend()
        axes[0].set_title("Radius Distribution")
        axes[0].set_xlabel("Radius"); axes[0].set_ylabel("Frequency")
        axes[0].grid(True, alpha=0.3)

        if Ref.shape[1] >= 3 and n_common > 0:
            idx = np.random.choice(n_common, min(n_common, 10_000), replace=False)
            sc  = axes[1].scatter(Ref[idx, 0], Ref[idx, 1],
                                  c=R_flat[idx], cmap="hot", s=1, alpha=0.8)
            plt.colorbar(sc, ax=axes[1], label="Radius")
            axes[1].set_aspect("equal"); axes[1].grid(True, alpha=0.3)
        axes[1].set_title("Reflector coloured by radius")
        axes[1].set_xlabel("X"); axes[1].set_ylabel("Y")

        plt.tight_layout()
        plt.savefig(j("reflector_analysis.png"), dpi=150, bbox_inches="tight")
        print("SAVED:", j("reflector_analysis.png"))

    # ── Density comparison ────────────────────────────────────────────────────
    if X_mesh is not None and Y_mesh is not None:
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))
        fig.suptitle("Source vs Destination Density", fontweight="bold")

        im1 = axes[0].imshow(X_mesh, extent=[-0.6, 0.6, -0.6, 0.6],
                             origin="lower", cmap="viridis")
        plt.colorbar(im1, ax=axes[0], label="Density")
        axes[0].set_title("Source"); axes[0].set_xlabel("X"); axes[0].set_ylabel("Y")

        im2 = axes[1].imshow(Y_mesh, extent=[-0.6, 0.6, -0.6, 0.6],
                             origin="lower", cmap="plasma")
        plt.colorbar(im2, ax=axes[1], label="Density")
        axes[1].set_title("Destination"); axes[1].set_xlabel("X"); axes[1].set_ylabel("Y")

        plt.tight_layout()
        plt.savefig(j("density_comparison.png"), dpi=150, bbox_inches="tight")
        print("SAVED:", j("density_comparison.png"))


if __name__ == "__main__":
    main()
