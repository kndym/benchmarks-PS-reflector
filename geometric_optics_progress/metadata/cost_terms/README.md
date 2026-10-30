# EOT terms and hard-map OT bounds

Updated October 30, 2026. Original kappa=0.6: 120 completed cases per map.
Repaired kappa=0.8: 72 completed cases per map (37 joint, zero fixed-epsilon,
35 fixed-source). The plots use the same frozen cases as the write-up's
raw EOT and debiased divergence figures. Pending cases are marked.

| Run | Map | EOT components | Hard-map OT bound |
| --- | --- | --- | --- |
| Original kappa=0.6 | Hard | [Plot](../../images/original_hard_cost_terms.png) | [Plot](../../images/original_hard_ot_bounds.png) |
| Original kappa=0.6 | Entropic + Snell | [Plot](../../images/original_entropic_cost_terms.png) | Not applicable |
| Repaired kappa=0.8 | Hard | [Plot](../../images/repaired_kappa_hard_cost_terms.png) | [Plot](../../images/repaired_kappa_hard_ot_bounds.png) |
| Repaired kappa=0.8 | Entropic + Snell | [Plot](../../images/repaired_kappa_entropic_cost_terms.png) | Not applicable |

Each component figure includes raw cross-EOT, push-forward self-EOT A, target
self-EOT B, their sum A+B, and the half-sum actually subtracted in Sinkhorn
divergence. This is the independent squared Euclidean EOT evaluation at
regularization 1/320 with KL relative to the product measure. The self terms
are distinct from the optical source/target corrections used to build the map.

For each hard map, pushed endpoints are matched exactly to target atoms and
source weights are aggregated into m. Keep s_i=min(m_i,q_i) at zero cost;
e_i=(m_i-q_i)+ and d_j=(q_j-m_j)+ leave common residual mass M. The requested
bound is U=sum_i e_i max_j ||y_i-y_j||^2. The figure also shows the tighter
maximum over deficit targets and the cost of the feasible coupling
pi_ij=s_i delta_ij+e_i d_j/M. All three bound unregularized OT; they are not
upper bounds on regularized EOT. Exact mass matching gives a zero OT bound.

All 192 hard cases passed witness marginal and bound-order checks; the largest
marginal discrepancy is 5.56e-17. An independent linear-programming test checks
that the bounds enclose the optimal unregularized transport cost on a small
nonuniform example. Other tests cover exact matching, empty target atoms,
concentrated pushed mass, and invalid weights.

[Manifest with source hashes](manifest.json). The six CSV files in this folder
contain the plotted values; the manifest identifies their hashes and saved-map
provenance. Regenerate with `python scripts/plot_writeup_cost_terms.py`.

[Write-up](../../README.md)
