"""
compare.py — Compare C++ (output_fast/) and Python (output_python/) results.

Run:  python compare.py

Loads:
  output_fast/R.txt, f.txt, g.txt, Ref.txt   (C++ plain-text output)
  output_python/R.npy, f.npy, g.npy, Ref.npy  (Python NumPy output)

Prints a numerical comparison report and shows side-by-side 3-D scatter plots.
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

REPO_ROOT   = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(REPO_ROOT, 'results')

# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

def load_txt_vec(path: str) -> np.ndarray:
    """Read a C++ output file: first line = N, then N scientific-notation values."""
    with open(path) as fh:
        lines = [l.strip() for l in fh if l.strip()]
    n = int(lines[0])
    vals = [float(v) for v in lines[1: n + 1]]
    if len(vals) < n:
        raise ValueError(f"{path}: expected {n} values, got {len(vals)}")
    return np.array(vals, dtype=np.float64)


def load_txt_mat(path: str, ncols: int = 3) -> np.ndarray:
    """Read a C++ matrix output: first line = N, then N lines of ncols values."""
    with open(path) as fh:
        lines = [l.strip() for l in fh if l.strip()]
    n = int(lines[0])
    rows = []
    for line in lines[1: n + 1]:
        rows.append([float(v) for v in line.split()])
    arr = np.array(rows, dtype=np.float64)
    if arr.shape != (n, ncols):
        raise ValueError(f"{path}: shape {arr.shape} != ({n}, {ncols})")
    return arr


# ---------------------------------------------------------------------------
# Load C++ results
# ---------------------------------------------------------------------------

cpp_dir = os.path.join(RESULTS_DIR, "output_fast")
py_dir  = os.path.join(RESULTS_DIR, "output_python")

missing = []
for d in [cpp_dir, py_dir]:
    if not os.path.isdir(d):
        missing.append(d)
if missing:
    print(f"ERROR: missing output directories: {missing}")
    print("Run  ./main_fast  and  python run_fast.py  first.")
    sys.exit(1)

print(f"Loading C++ results from {cpp_dir}/")
try:
    R_cpp   = load_txt_vec(f"{cpp_dir}/R.txt")
    f_cpp   = load_txt_vec(f"{cpp_dir}/f.txt")
    g_cpp   = load_txt_vec(f"{cpp_dir}/g.txt")
    Ref_cpp = load_txt_mat(f"{cpp_dir}/Ref.txt", ncols=3)
    f_id_cpp = load_txt_vec(f"{cpp_dir}/f_id.txt")
    g_id_cpp = load_txt_vec(f"{cpp_dir}/g_id.txt")
except Exception as e:
    print(f"ERROR reading C++ output: {e}")
    sys.exit(1)

print(f"Loading Python results from {py_dir}/")
try:
    R_py    = np.load(f"{py_dir}/R.npy")
    f_py    = np.load(f"{py_dir}/f.npy")
    g_py    = np.load(f"{py_dir}/g.npy")
    Ref_py  = np.load(f"{py_dir}/Ref.npy")
    f_id_py = np.load(f"{py_dir}/f_id.npy")
    g_id_py = np.load(f"{py_dir}/g_id.npy")
except Exception as e:
    print(f"ERROR reading Python output: {e}")
    sys.exit(1)

# ---------------------------------------------------------------------------
# Numerical comparison
# ---------------------------------------------------------------------------

NK_cpp = len(R_cpp)
NK_py  = len(R_py)

print("\n" + "=" * 60)
print(f"  NK (C++) = {NK_cpp},  NK (Python) = {NK_py}")
if NK_cpp != NK_py:
    print("  WARNING: different sizes — comparison may be invalid!")
    sys.exit(1)
NK = NK_cpp
print("=" * 60)

def report(name, a, b):
    diff     = np.abs(a - b)
    rel      = diff / (np.abs(a) + 1e-10)
    print(f"\n  {name}:")
    print(f"    max |diff|   = {diff.max():.6e}")
    print(f"    mean |diff|  = {diff.mean():.6e}")
    print(f"    max rel err  = {rel.max():.6e}")
    return diff.max()

print("\nArray-wise comparison:")
max_R   = report("R   (reflector radius)", R_cpp,   R_py)
max_f   = report("f   (source potential)", f_cpp,   f_py)
max_g   = report("g   (target potential)", g_cpp,   g_py)
report(           "f_id",                   f_id_cpp, f_id_py)
report(           "g_id",                   g_id_cpp, g_id_py)

print("\n" + "-" * 60)
print(f"  R summary:")
print(f"    C++:    min={R_cpp.min():.4f}, max={R_cpp.max():.4f}, mean={R_cpp.mean():.4f}")
print(f"    Python: min={R_py.min():.4f},  max={R_py.max():.4f},  mean={R_py.mean():.4f}")
print("-" * 60)

if max_R < 0.05:
    print(f"\n  PASS: max|R_cpp - R_py| = {max_R:.4e} < 0.05")
else:
    print(f"\n  WARN: max|R_cpp - R_py| = {max_R:.4e} >= 0.05  — larger than expected")

# ---------------------------------------------------------------------------
# Side-by-side 3-D scatter plots
# ---------------------------------------------------------------------------

try:
    from mpl_toolkits.mplot3d import Axes3D   # noqa: F401

    vmin = min(R_cpp.min(), R_py.min())
    vmax = max(R_cpp.max(), R_py.max())

    fig = plt.figure(figsize=(14, 6))

    ax1 = fig.add_subplot(121, projection="3d")
    sc1 = ax1.scatter(Ref_cpp[:, 0], Ref_cpp[:, 1], Ref_cpp[:, 2],
                      c=R_cpp, cmap="viridis", s=15, vmin=vmin, vmax=vmax)
    plt.colorbar(sc1, ax=ax1, label="R")
    ax1.set_title(f"C++  (NK={NK})")

    ax2 = fig.add_subplot(122, projection="3d")
    sc2 = ax2.scatter(Ref_py[:, 0], Ref_py[:, 1], Ref_py[:, 2],
                      c=R_py,  cmap="viridis", s=15, vmin=vmin, vmax=vmax)
    plt.colorbar(sc2, ax=ax2, label="R")
    ax2.set_title(f"Python  (NK={NK})")

    fig.suptitle("Reflector surface comparison", fontsize=13)
    plt.tight_layout()
    plt.show()

except Exception as e:
    print(f"[plot skipped: {e}]")
