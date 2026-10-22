# Exploratory density research

All computations import the active repository pipeline; no repository solver changes. Uniform source and existing uniform, gaussian, donut and cross target factories use the original source/target spherical patches and kappa=0.6. Halton samples are nested. The main sweep uses original warm-start run_sinkhorn_divergence and hard c-transform pushforwards. PF evaluation uses squared Euclidean cost in 3D.

Completed 192 PF evaluations. See density_costs.csv; finite original cost does not certify marginal feasibility.

## cross
Original solver untruncated marginal L1 residual range: 0.00147036 to 0.0177877.
Refractor eps=0.000195313; NK in increasing-cost order: PF eps=0.000195313: 400,225,100; PF eps=0.00078125: 400,225,100; PF eps=0.003125: 400,225,100; PF eps=0.0125: 400,225,100.
Refractor eps=0.00078125; NK in increasing-cost order: PF eps=0.000195313: 225,400,100; PF eps=0.00078125: 225,400,100; PF eps=0.003125: 225,100,400; PF eps=0.0125: 225,100,400.
Refractor eps=0.003125; NK in increasing-cost order: PF eps=0.000195313: 100,225,400; PF eps=0.00078125: 100,225,400; PF eps=0.003125: 100,225,400; PF eps=0.0125: 100,225,400.
Refractor eps=0.0125; NK in increasing-cost order: PF eps=0.000195313: 225,400,100; PF eps=0.00078125: 225,400,100; PF eps=0.003125: 225,100,400; PF eps=0.0125: 225,100,400.

## donut
Original solver untruncated marginal L1 residual range: 0.00335298 to 0.0277039.
Refractor eps=0.000195313; NK in increasing-cost order: PF eps=0.000195313: 400,225,100; PF eps=0.00078125: 400,225,100; PF eps=0.003125: 400,225,100; PF eps=0.0125: 225,400,100.
Refractor eps=0.00078125; NK in increasing-cost order: PF eps=0.000195313: 225,400,100; PF eps=0.00078125: 100,225,400; PF eps=0.003125: 225,400,100; PF eps=0.0125: 225,400,100.
Refractor eps=0.003125; NK in increasing-cost order: PF eps=0.000195313: 100,225,400; PF eps=0.00078125: 100,225,400; PF eps=0.003125: 100,225,400; PF eps=0.0125: 100,225,400.
Refractor eps=0.0125; NK in increasing-cost order: PF eps=0.000195313: 225,400,100; PF eps=0.00078125: 225,400,100; PF eps=0.003125: 225,400,100; PF eps=0.0125: 225,100,400.

## gaussian
Original solver untruncated marginal L1 residual range: 0.00128318 to 0.0181276.
Refractor eps=0.000195313; NK in increasing-cost order: PF eps=0.000195313: 400,225,100; PF eps=0.00078125: 400,225,100; PF eps=0.003125: 400,225,100; PF eps=0.0125: 400,225,100.
Refractor eps=0.00078125; NK in increasing-cost order: PF eps=0.000195313: 225,400,100; PF eps=0.00078125: 225,400,100; PF eps=0.003125: 225,400,100; PF eps=0.0125: 400,225,100.
Refractor eps=0.003125; NK in increasing-cost order: PF eps=0.000195313: 100,225,400; PF eps=0.00078125: 100,225,400; PF eps=0.003125: 100,225,400; PF eps=0.0125: 100,225,400.
Refractor eps=0.0125; NK in increasing-cost order: PF eps=0.000195313: 225,400,100; PF eps=0.00078125: 225,400,100; PF eps=0.003125: 225,100,400; PF eps=0.0125: 225,100,400.

## uniform
Original solver untruncated marginal L1 residual range: 0.00382488 to 0.0305634.
Refractor eps=0.000195313; NK in increasing-cost order: PF eps=0.000195313: 100,225,400; PF eps=0.00078125: 100,400,225; PF eps=0.003125: 100,400,225; PF eps=0.0125: 400,225,100.
Refractor eps=0.00078125; NK in increasing-cost order: PF eps=0.000195313: 225,100,400; PF eps=0.00078125: 100,225,400; PF eps=0.003125: 100,225,400; PF eps=0.0125: 100,225,400.
Refractor eps=0.003125; NK in increasing-cost order: PF eps=0.000195313: 225,100,400; PF eps=0.00078125: 100,225,400; PF eps=0.003125: 100,225,400; PF eps=0.0125: 100,225,400.
Refractor eps=0.0125; NK in increasing-cost order: PF eps=0.000195313: 225,100,400; PF eps=0.00078125: 100,225,400; PF eps=0.003125: 100,225,400; PF eps=0.0125: 100,225,400.

## Convergence control and limitations
A separate dense untruncated Sinkhorn refinement starts from original potentials and alternates exact logsumexp updates until max(source L1 error,target L1 error)<1e-8, capped at3000 iterations. This is research-only code in refine.py, not a proposed pipeline edit. refinement.csv compares maps/costs.

Concrete ranking reversal: at refractor epsilon=1/80 and PF epsilon=1/5120, original uniform costs are 0.0018612772 (NK100) and0.0018829629 (NK400). Refined costs are0.0019109488 (NK100) and0.0018332828 (NK400), reversing the ordering. Thus apparent NK order may reflect the inherited iteration/truncation schedule as well as discretization.

cold_probe.py/cold_probe.csv preserve a failed cold-start diagnostic: zero potentials with C++ threshold log(1e-6/NK) discard essentially all cost entries for these patches before updating, with full-plan L1 residual near1. Increasing cap from16 to512 does not fix this initialization. The cold probe is not interpreted as solved transport.

For a fixed PF and target measure, entropic objective OT_epsilon=<c,pi>+epsilon KL(pi|p tensor q) is nondecreasing in epsilon because KL>=0; this says nothing about NK ordering or ordering when PF changes. The hard argmin map can jump with potentials; log-domain numeric stability alone does not guarantee solved marginals. Epsilon below1/5120 and NK beyond400 were omitted to keep this side task bounded.

## Recorded objectives and certification
The source_objective field is the raw cross dual objective reconstructed from f_raw/g_raw, with the exponential mass correction. It is not the original total_cost, which subtracts approximate source/target identity potential terms. At infeasible marginals, source_objective is a dual bound and not a certified primal transport cost. All main solver runs retain inherited C++ cap/truncation semantics. PF costs are separately evaluated by entropic_ot_cost; PF convergence flags certify only that second computation under its current potential-change stopping rule.
PF evaluations reporting nonconvergence: 0. Source problems with marginal L1 residual above 1e-6: 48. Exact source certification uses independent residual checks, not the original code stopping report.