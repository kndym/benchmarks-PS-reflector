# Sinkhorn convergence research

This side investigation imports the active refractor solver and continues its existing `sinkhorn_step` outside the production pipeline. It changes no solver defaults or repository code. `probe.py` diagnoses 3 point counts (100, 225, 400), 4 epsilon values (0.0125 to 0.000125, ratio 100), and uniform versus smoothly tilted source/target densities on the original patches. Results are checkpointed in `diagnostics.json`; final results also appear in `diagnostics.csv`.

## Recommendation

For normalized balanced probability measures use

`residual = max(sum(abs(row_marginal - p)), sum(abs(column_marginal - q)))`

with **1e-6 as a practical default**, and **1e-8 for comparative research** where cost differences or rankings are small. Those numbers are recommendations, not universal optimal thresholds. Verify sensitive rankings by repeating at a tighter tolerance. Check every 5–10 iterations and always at final return. Log-domain marginal evaluation avoids constructing the full transport plan and can be chunked with the existing cost blocks. The marginals must use the raw cross-OT potentials and full cost, not identity-subtracted surface potentials, and should include kernel entries discarded by the iteration's screening rule.

For this API the implied plan is

`log(pi_ij) = log(p_i) + log(q_j) + (f_i + g_j - C_ij)/epsilon`.

Evaluate row/column sums with logsumexp, exponentiating only their final values. Log-domain arithmetic prevents underflow of the kernel from being mistaken for convergence. Zero-weight support must be excluded appropriately.

Keep a maximum iteration count, e.g. 2000, as a safety budget and return an explicit convergence status plus residual. A budget exit is not convergence. Particularly small epsilon can require more than 2000 iterations: it is better to report an incomplete solve than silently assign a cost ranking. Use warm starts and epsilon continuation, with more inner iterations where needed, instead of assuming one iteration per epsilon stage is sufficient.

For highly unequal weights also report the worst relative marginal error on sufficiently large positive atoms. L1 is stable as NK changes; a fixed L2 threshold can become easier to satisfy merely because the same mass discrepancy is spread across more atoms. Relative checks alone can overemphasize negligible atoms, so they should be an additional diagnostic.

Potential changes and objective stagnation are useful secondary diagnostics but should not replace feasibility. Potential gauges can change without changing the plan, and the potential change relevant to exponential scalings is change/epsilon. If cost accuracy is the real target, a feasible-primal versus dual gap is a stronger certificate; this requires rounding the plan to exact marginals and computing the entropy term, and is more work than a minimal default. A primal computed from the unbalanced current plan does not supply a valid upper bound.

## What is inherited here

The active compatibility loop allows 17 final iterations because it tests `i > cap_iter` with `cap_iter=16`. It checks `max(abs(log(F)/k))` against `1e-5`, after computing G then F, and absorbs half of each scaling update. Thus this diagnostic is epsilon times a source-side log scaling correction, evaluated before damped absorption. It does not check both actual plan marginals.

Keeping 1e-5 fixed means cap_thr/epsilon is 0.0008 at epsilon=0.0125 but 0.08 at epsilon=0.000125. As a rough scaling intuition those correspond to factors exp(0.0008) and exp(0.08). This is not a feasibility certificate because the checked F is based on an intermediate G and the absorption is damped.

The `entropic_ot_cost` PF evaluator is more rigorous about its objective (it includes the exponential dual term), but its stop uses absolute potential change 1e-9. That too should eventually gain an actual marginal residual criterion. No edits were made here.

The inherited identity routines use `C(x,y)` with source/source or target/target weights on the two indices. They do not evaluate conventional `C(x,x)` and `C(y,y)` self transports. Because x and y are different patches and C(x,y) is generally not symmetric as an indexed matrix, the symmetric single-potential iteration cannot be assumed to solve standard self OT. Preserve these corrections for C++ compatibility, but avoid describing their aggregate as a certified conventional Sinkhorn divergence. This investigation measures convergence of the raw cross plan independently.

## Evidence and scope

All 24 planned cases completed. The maximum inherited L1 residual was 0.03056. Six of 24 cases remained above 1e-6 after the extra 2000-iteration budget, all at epsilon=0.000125. The summary figure is `inherited_residual.png`.

At NK=225 with uniform densities, the inherited result had L1 residual 0.01813 for epsilon=0.0125 and 0.01123 for epsilon=0.003125. Continuing the exact imported step needed another 40 and 140 iterations, respectively, to reach 1e-6. At epsilon=0.00078125 the inherited residual was 0.00715 and another 500 iterations were required. The inherited potential threshold first became satisfied while residuals remained 0.000421, 0.001998, and 0.007027, respectively. These are observed results, not inferred error bounds.

At NK=100, tilted weights, epsilon=0.000125, even an extra 2000 inherited iterations left residual 4.12e-6. All numerical quantities stayed finite. Finite arithmetic and near-unit total mass therefore do not establish accurate marginals. The dual cost moved only 4.08e-7 while feasibility was visibly worse, illustrating why objective stability alone is insufficient.

## Primary sources

- [POT's actual Sinkhorn implementation](https://pythonot.github.io/_modules/ot/bregman/_sinkhorn.html) evaluates marginal violations periodically (every 10 iterations in its log-domain implementation) and reports iteration exhaustion separately. Its epsilon-scaling implementation checks both marginals. POT uses a norm and tolerance different from the explicit L1 recommendation above.
- [Schmitzer, Stabilized Sparse Scaling Algorithms for Entropy Regularized Transport Problems](https://arxiv.org/abs/1610.06519) motivates stabilization, epsilon scaling, and sparsity as small epsilon slows convergence and makes direct kernel arithmetic unstable.
- [Altschuler, Weed, Rigollet, Near-linear time approximation algorithms for optimal transport via Sinkhorn iteration](https://arxiv.org/abs/1705.09634) provides a theoretical framework separating marginal approximation from feasible transport construction.
- The local Benamou–Ijzerman–Rukhaia paper, equations (2.1), (2.2), and (2.5)–(2.6), defines entropic duals and genuine diagonal transports. Extracted text is retained as `paper.txt` solely for this local investigation.
