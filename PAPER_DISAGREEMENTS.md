# Python/C++ pipeline differences intentionally deferred

These items were found during the Python-to-C++ and paper audit, but are not
changed in this pass because they require a separate decision about whether
the paper or the original C++ implementation is authoritative.

## Diagonal Sinkhorn-divergence problems

Paper equations (2.5)--(2.7) define the debiasing terms using the diagonal
problems `(mu, mu)` and `(nu, nu)`.  That requires the source identity solve to
use `C(x_i, x_j)` and the target identity solve to use `C(y_i, y_j)`.

The original C++ code and the Python translation instead use the cross-grid
cost `C(x_i, y_j)` for both identity solves.  This behavior remains unchanged
for C++ parity.  It is a paper-level algorithm discrepancy that affects the
meaning of the reported Sinkhorn divergence.

## Sinkhorn update ordering

Paper equation (2.2) presents the block update as `f` first and then `g`.
The original C++ `Sinkhorn_axb` routine updates `G` first and then `F`.
Python follows the C++ order.  Both orders have the same fixed point, but a
finite, capped run can produce different intermediate potentials.

## Entropic interpolation

The paper's equations (3.7)--(3.8) define soft-min entropic interpolations
for the OT and Sinkhorn-divergence potentials.  The current pipeline uses the
hard minimum c-transform (paper equation (3.6)) for the regular-grid
reflector.  Adding the soft-min interpolation is deferred.

## Paper ray-tracing methodology

The paper evaluates the continuous interpolated reflector by ray tracing.  A
few convenience scripts in this repository use a discrete argmin OT map for
their pushed points instead.  The C++ reflector executable uses the regular
grid ray tracer; replacing the convenience-script map with the full paper
interpolation/ray-tracing workflow is deferred.
