# Uniform refractor sensitivity controls

Ten saved controls vary kappa or translate one source-patch angle while keeping
the remaining baseline parameters fixed. All use 1600 source/target samples.
The exact parameters and current diagnostics are in
[`L2_costs_numeric.csv`](L2_costs_numeric.csv),
[`metadata.json`](metadata.json), and [`L2_costs.txt`](L2_costs.txt).

## Current metric definitions

The current summary reports debiased entropic Sinkhorn divergence between the
pushed unit-sphere points and the prescribed target, with **squared Euclidean
3D ground cost and epsilon 0.005**. It also reports its square root and the
convergence/iteration diagnostics for cross and self transports. This is not
the non-debiased EOT at epsilon 1/320 used in the newer parameter sweeps, and
its square root is not exact unregularized Wasserstein distance.

The separate mass-L2 diagnostic compares pushed and prescribed weights at the
target atoms. A passed PF stopping test does not certify marginal feasibility
of the earlier optical source solve.

## Surface images

The laptop's `*_surface.png` images render the corresponding existing case NPZ
surface arrays; they add visualizations, not ten new numerical experiments.
The filenames identify their case. The surface image manifest records the
source import commit and hashes linking each image to its local NPZ.

The NPZs and current CSV/metadata already existed on the desktop branch. They
are retained unchanged when importing these laptop-only renders. Earlier
laptop reports may describe historical metrics; use the current CSV/metadata
for the definitions above.

## Case gallery

| Case | Kappa | Surface image | Existing numerical result |
| --- | ---: | --- | --- |
| original_uniform | 0.6 | [View](original_uniform_surface.png) | [NPZ](original_uniform.npz) |
| kappa_0p25 | 0.25 | [View](kappa_0p25_surface.png) | [NPZ](kappa_0p25.npz) |
| kappa_0p50 | 0.5 | [View](kappa_0p50_surface.png) | [NPZ](kappa_0p50.npz) |
| kappa_0p90 | 0.9 | [View](kappa_0p90_surface.png) | [NPZ](kappa_0p90.npz) |
| source_phi_plus_05deg | 0.6 | [View](source_phi_plus_05deg_surface.png) | [NPZ](source_phi_plus_05deg.npz) |
| source_phi_plus_15deg | 0.6 | [View](source_phi_plus_15deg_surface.png) | [NPZ](source_phi_plus_15deg.npz) |
| source_phi_plus_25deg | 0.6 | [View](source_phi_plus_25deg_surface.png) | [NPZ](source_phi_plus_25deg.npz) |
| source_theta_plus_10deg | 0.6 | [View](source_theta_plus_10deg_surface.png) | [NPZ](source_theta_plus_10deg.npz) |
| source_theta_plus_25deg | 0.6 | [View](source_theta_plus_25deg_surface.png) | [NPZ](source_theta_plus_25deg.npz) |
| source_theta_plus_40deg | 0.6 | [View](source_theta_plus_40deg_surface.png) | [NPZ](source_theta_plus_40deg.npz) |

Import and hash checks: [desktop validation](../archive/laptop_research/desktop_import_validation.json).
