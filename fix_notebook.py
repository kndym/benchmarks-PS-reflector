"""
fix_notebook.py  —  patches the two bugs in Benchmark.ipynb's VIZ_SCRIPT
"""
import json, os

NB_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Benchmark.ipynb')

with open(NB_PATH, 'r', encoding='utf-8') as f:
    nb = json.load(f)

for cell in nb['cells']:
    if cell.get('id') == 'viz-cell':
        src = ''.join(cell['source'])

        # ── Fix 1: numpy array `or` fallback ─────────────────────────────────
        OLD_OR = 'Y_push = load_points(j("Y_Pushed_projected.txt")) or load_points(j("Y_Pushed_projected_owndisc.txt"))\n'
        NEW_OR = (
            '_yp = load_points(j("Y_Pushed_projected.txt"))\n'
            'Y_push = _yp if _yp is not None else load_points(j("Y_Pushed_projected_owndisc.txt"))\n'
        )
        if OLD_OR in src:
            src = src.replace(OLD_OR, NEW_OR, 1)
            print("OK: Fixed Y_push numpy `or` -> explicit None check")
        else:
            print("WARN: Fix 1 pattern not found")
            print("  Near Y_push:", repr(src[src.find('Y_push'):src.find('Y_push')+120]))

        # ── Fix 2: Axes3D import guard ────────────────────────────────────────
        OLD_3D = (
            'ax = plt.subplot(2, 3, 5, projection="3d")\n'
            'if Ref is not None and len(Ref) and Ref.shape[1] >= 3:\n'
            '    idx = np.random.choice(len(Ref), min(len(Ref), 10000), replace=False)\n'
            '    ax.scatter(Ref[idx,0], Ref[idx,1], Ref[idx,2], s=0.5, alpha=0.6, c=Ref[idx,2], cmap="viridis")\n'
            '    ax.set_title("Reflector Surface", fontweight="bold")\n'
            '    ax.set_xlabel("X"); ax.set_ylabel("Y"); ax.set_zlabel("Z")\n'
            'else: ax.text2D(0.5, 0.5, "Ref_MY.txt not found", ha="center", va="center", transform=ax.transAxes)\n'
        )
        NEW_3D = (
            'try:\n'
            '    from mpl_toolkits.mplot3d import Axes3D\n'
            '    _has_3d = True\n'
            'except Exception:\n'
            '    _has_3d = False\n'
            'if _has_3d:\n'
            '    ax = plt.subplot(2, 3, 5, projection="3d")\n'
            '    if Ref is not None and len(Ref) and Ref.shape[1] >= 3:\n'
            '        idx = np.random.choice(len(Ref), min(len(Ref), 10000), replace=False)\n'
            '        ax.scatter(Ref[idx,0], Ref[idx,1], Ref[idx,2], s=0.5, alpha=0.6, c=Ref[idx,2], cmap="viridis")\n'
            '        ax.set_title("Reflector Surface (3D)", fontweight="bold")\n'
            '        ax.set_xlabel("X"); ax.set_ylabel("Y"); ax.set_zlabel("Z")\n'
            '    else:\n'
            '        ax.text2D(0.5, 0.5, "Ref_MY.txt not found", ha="center", va="center", transform=ax.transAxes)\n'
            'else:\n'
            '    ax = plt.subplot(2, 3, 5)\n'
            '    if Ref is not None and len(Ref) and Ref.shape[1] >= 3:\n'
            '        idx = np.random.choice(len(Ref), min(len(Ref), 10000), replace=False)\n'
            '        sc = ax.scatter(Ref[idx,0], Ref[idx,1], c=Ref[idx,2], cmap="viridis", s=0.5, alpha=0.6)\n'
            '        plt.colorbar(sc, ax=ax, label="Z")\n'
            '        ax.set_title("Reflector Surface (2D, Z=colour)", fontweight="bold")\n'
            '        ax.set_aspect("equal"); ax.grid(True, alpha=0.3)\n'
            '    else:\n'
            '        no_data(ax, "Ref_MY.txt not found")\n'
            'ax.set_xlabel("X"); ax.set_ylabel("Y")\n'
        )

        # The shape[1] comparison may be stored as \u003e= in JSON — normalise
        OLD_3D_ALT = OLD_3D.replace('>=', '\u003e=')
        NEW_3D_ALT = NEW_3D.replace('>=', '\u003e=')

        if OLD_3D in src:
            src = src.replace(OLD_3D, NEW_3D, 1)
            print("OK: Fixed Axes3D -> try/except with 2D fallback")
        elif OLD_3D_ALT in src:
            src = src.replace(OLD_3D_ALT, NEW_3D_ALT, 1)
            print("OK: Fixed Axes3D (unicode variant) -> try/except with 2D fallback")
        else:
            print("WARN: Fix 2 pattern not found")
            i3d = src.find('subplot(2, 3, 5')
            print("  Context:", repr(src[i3d:i3d+200]))

        cell['source'] = [src]
        break
else:
    raise RuntimeError("viz-cell not found in notebook")

with open(NB_PATH, 'w', encoding='utf-8') as f:
    json.dump(nb, f, indent=1, ensure_ascii=False)

print("\nNotebook saved:", NB_PATH)
