# Compare C++ and Python reflector outputs for one requested point-cloud size.
# Usage: python compare_nk.py 300

import os, sys
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(SCRIPT_DIR))
FIGURES_DIR = os.path.join(REPO_ROOT, "figures")
os.makedirs(FIGURES_DIR, exist_ok=True)

# Select the output directories for the requested point-cloud size.
NK_arg = int(sys.argv[1]) if len(sys.argv) > 1 else 300
cpp_dir = f"output_cpp_NK{NK_arg}"
py_dir  = f"output_py_NK{NK_arg}"

for d in [cpp_dir, py_dir]:
    if not os.path.isdir(d):
        print(f"ERROR: directory not found: {d}")
        print("Run  main_compare_{NK}  and  python run_compare.py {NK}  first.")
        sys.exit(1)

# Read the vector and surface arrays used in the comparison.
def load_vec(path):
    with open(path) as fh:
        lines = [l.strip() for l in fh if l.strip()]
    n = int(lines[0])
    return np.array([float(v) for v in lines[1:n+1]], dtype=np.float64)

def load_mat(path, ncols=3):
    with open(path) as fh:
        lines = [l.strip() for l in fh if l.strip()]
    n = int(lines[0])
    rows = [[float(v) for v in line.split()] for line in lines[1:n+1]]
    return np.array(rows, dtype=np.float64)

print(f"Loading C++ from {cpp_dir}/")
R_cpp    = load_vec(f"{cpp_dir}/R.txt")
f_cpp    = load_vec(f"{cpp_dir}/f.txt")
g_cpp    = load_vec(f"{cpp_dir}/g.txt")
f_id_cpp = load_vec(f"{cpp_dir}/f_id.txt")
g_id_cpp = load_vec(f"{cpp_dir}/g_id.txt")
Ref_cpp  = load_mat(f"{cpp_dir}/Ref.txt")

print(f"Loading Python from {py_dir}/")
R_py     = np.load(f"{py_dir}/R.npy")
f_py     = np.load(f"{py_dir}/f.npy")
g_py     = np.load(f"{py_dir}/g.npy")
f_id_py  = np.load(f"{py_dir}/f_id.npy")
g_id_py  = np.load(f"{py_dir}/g_id.npy")
Ref_py   = np.load(f"{py_dir}/Ref.npy")

NK_cpp, NK_py = len(R_cpp), len(R_py)
print(f"\n{'='*58}")
print(f"  NK (C++) = {NK_cpp},  NK (Python) = {NK_py}")
if NK_cpp != NK_py:
    print("  WARNING: sizes differ — aborting comparison.")
    sys.exit(1)
NK = NK_cpp
print(f"{'='*58}")

# Report absolute and relative differences between matching arrays.
def report(name, a, b):
    diff = np.abs(a - b)
    rel  = diff / (np.abs(a) + 1e-10)
    print(f"\n  {name}:")
    print(f"    max |diff|  = {diff.max():.6e}")
    print(f"    mean |diff| = {diff.mean():.6e}")
    print(f"    max rel err = {rel.max():.6e}")
    return diff.max()

print("\nArray-wise comparison:")
max_R    = report("R   (reflector radius)", R_cpp,    R_py)
report(           "f   (source potential)", f_cpp,    f_py)
report(           "g   (target potential)", g_cpp,    g_py)
report(           "f_id",                  f_id_cpp,  f_id_py)
report(           "g_id",                  g_id_cpp,  g_id_py)

print(f"\n{'-'*58}")
print(f"  R (C++):   min={R_cpp.min():.4f}  max={R_cpp.max():.4f}  mean={R_cpp.mean():.4f}")
print(f"  R (Python):min={R_py.min():.4f}  max={R_py.max():.4f}  mean={R_py.mean():.4f}")
print(f"{'-'*58}")

thresh = 0.05
if max_R < thresh:
    print(f"\n  PASS: max|R_cpp - R_py| = {max_R:.4e} < {thresh}")
else:
    print(f"\n  WARN: max|R_cpp - R_py| = {max_R:.4e} >= {thresh}")

# Save a side-by-side 3-D surface comparison.
try:
    from mpl_toolkits.mplot3d import Axes3D  # noqa
    vmin = min(R_cpp.min(), R_py.min())
    vmax = max(R_cpp.max(), R_py.max())

    fig = plt.figure(figsize=(14, 6))
    ax1 = fig.add_subplot(121, projection="3d")
    sc1 = ax1.scatter(Ref_cpp[:,0], Ref_cpp[:,1], Ref_cpp[:,2],
                      c=R_cpp, cmap="viridis", s=15, vmin=vmin, vmax=vmax)
    plt.colorbar(sc1, ax=ax1, label="R")
    ax1.set_title(f"C++  NK={NK}")

    ax2 = fig.add_subplot(122, projection="3d")
    sc2 = ax2.scatter(Ref_py[:,0], Ref_py[:,1], Ref_py[:,2],
                      c=R_py, cmap="viridis", s=15, vmin=vmin, vmax=vmax)
    plt.colorbar(sc2, ax=ax2, label="R")
    ax2.set_title(f"Python  NK={NK}")

    fig.suptitle(f"Reflector surface comparison — NK={NK}", fontsize=13)
    plt.tight_layout()
    fig_path = os.path.join(FIGURES_DIR, f"fig_compare_NK{NK}.png")
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.show()
    print(f"\nSaved {fig_path}")
except Exception as e:
    print(f"[plot skipped: {e}]")
