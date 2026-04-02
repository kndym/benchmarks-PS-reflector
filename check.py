



import subprocess, importlib

try:
    import plotly.graph_objects as go
except ImportError:
    subprocess.run([sys.executable, '-m', 'pip', 'install', 'plotly', '--quiet'], check=True)
    importlib.invalidate_caches()
    import plotly.graph_objects as go

MAX_PTS = 6_000
rng     = np.random.default_rng(42)

fig3d = go.Figure()

# ── 1. Source hemisphere: scatter coloured by P(x) ──────────────────────────
src_idx = np.where(mask_p)[0]
si      = rng.choice(len(src_idx), min(len(src_idx), MAX_PTS), replace=False)
s_pts   = x[src_idx[si]];  s_col = p[src_idx[si]]
fig3d.add_trace(go.Scatter3d(
    x=s_pts[:, 0], y=s_pts[:, 1], z=s_pts[:, 2],
    mode='markers',
    marker=dict(size=2.5, color=s_col, colorscale='YlOrRd', opacity=0.65, showscale=False),
    name='Source  P(x)',
    hovertemplate='x=%{x:.3f}  y=%{y:.3f}  z=%{z:.3f}<br>P=%{marker.color:.4f}<extra></extra>',
))

# ── 2. Refractor surface: scatter coloured by radius ─────────────────────────
ref_idx = np.where(mask_p)[0]
ri      = rng.choice(len(ref_idx), min(len(ref_idx), MAX_PTS), replace=False)
ref_pts = Ref[ref_idx[ri]];  r_col = R[ref_idx[ri]]
fig3d.add_trace(go.Scatter3d(
    x=ref_pts[:, 0], y=ref_pts[:, 1], z=ref_pts[:, 2],
    mode='markers',
    marker=dict(size=3, color=r_col, colorscale='Viridis', opacity=0.85,
                showscale=True,
                colorbar=dict(title='Radius', x=1.02, len=0.55, y=0.5, thickness=15)),
    name='Refractor surface',
    hovertemplate='x=%{x:.3f}  y=%{y:.3f}  z=%{z:.3f}<br>R=%{marker.color:.4f}<extra></extra>',
))

# ── 3. Output light: push-forward on upper hemisphere ───────────────────────
n_push = len(y_push_3d)
pi     = rng.choice(n_push, min(n_push, MAX_PTS), replace=False)
push_col = Y_Pushed_projected[pi, 2]   # Q(y_pushed)
fig3d.add_trace(go.Scatter3d(
    x=y_push_3d[pi, 0], y=y_push_3d[pi, 1], z=y_push_3d[pi, 2],
    mode='markers',
    marker=dict(size=2.5, color=push_col, colorscale='Blues', opacity=0.65, showscale=False),
    name='Output light  Q(y)',
    hovertemplate='x=%{x:.3f}  y=%{y:.3f}  z=%{z:.3f}<br>Q=%{marker.color:.4f}<extra></extra>',
))

fig3d.update_layout(
    title=dict(
        text=(
            '<b>Optical System -- Refraction</b><br>'
            f'<sup>NK={NK}  \u00b7  \u03ba={kappa}  \u00b7  '
            'source P(x)  \u00b7  refractor (radius)  \u00b7  output Q(y)  \u00b7  '
            'drag=rotate  scroll=zoom  hover=values</sup>'
        ),
        x=0.5, xanchor='center',
    ),
    scene=dict(
        xaxis=dict(title='X', showgrid=True),
        yaxis=dict(title='Y', showgrid=True),
        zaxis=dict(title='Z', showgrid=True),
        aspectmode='data',
        camera=dict(eye=dict(x=1.4, y=1.4, z=0.8)),
    ),
    legend=dict(itemsizing='constant', x=0, y=1),
    margin=dict(l=0, r=80, t=100, b=0),
    width=960, height=720,
)

html_path = 'fig_optical_3d.html'
fig3d.write_html(html_path, include_plotlyjs='cdn')
from IPython.display import HTML, display as _display
_display(HTML(filename=html_path))
print(f'Interactive 3D plot saved to {html_path}')