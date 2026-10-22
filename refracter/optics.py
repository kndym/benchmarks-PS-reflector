"""Normals of the entropic radial surface and vector Snell refraction."""
import numpy as np
from scipy.special import logsumexp


def _unit_vectors(value, name):
    value = np.asarray(value, dtype=float)
    if value.ndim != 2 or value.shape[1] != 3 or not np.isfinite(value).all():
        raise ValueError(f"{name} must be finite (N, 3) vectors")
    if not np.allclose(np.linalg.norm(value, axis=1), 1, atol=1e-10, rtol=0):
        raise ValueError(f"{name} must contain unit vectors")
    return value


def entropic_surface_normals(x, y, q, g, epsilon, kappa, chunk_size=512):
    """Return unnormalized normals N and log R for R=exp(softmin(g-c)).

    N=sum_j w_j (x-kappa*y_j)/(1-kappa*x.y_j), with normalized softmin
    weights q_j exp((g_j-c_j)/epsilon). Uses raw target potential g,
    not the self-corrected potential. N satisfies x.N=1. The log radius
    avoids overflow from the arbitrary additive gauge of g.
    """
    x, y = _unit_vectors(x, "x"), _unit_vectors(y, "y")
    q, g = np.asarray(q, dtype=float), np.asarray(g, dtype=float)
    if q.shape != (len(y),) or g.shape != q.shape:
        raise ValueError("q and g must match the target point count")
    if not np.isfinite(q).all() or (q < 0).any() or not q.sum() > 0:
        raise ValueError("q must be finite nonnegative weights with positive sum")
    if not np.isfinite(g[q > 0]).all():
        raise ValueError("Supported target potentials must be finite")
    if not np.isfinite(epsilon) or epsilon <= 0 or not 0 < kappa < 1 or chunk_size < 1:
        raise ValueError("Require epsilon>0, 0<kappa<1, and positive chunk size")
    y, g, q = y[q > 0], g[q > 0], q[q > 0] / q.sum()
    normals, log_radius = np.empty_like(x), np.empty(len(x))
    for start in range(0, len(x), chunk_size):
        xx = x[start:start + chunk_size]
        denominator = 1 - kappa * (xx @ y.T)
        logits = np.log(q)[None, :] + (g[None, :] + np.log(denominator)) / epsilon
        log_z = logsumexp(logits, axis=1)
        weighted = np.exp(logits - log_z[:, None]) / denominator
        normals[start:start + len(xx)] = (xx * weighted.sum(axis=1)[:, None]
                                          - kappa * (weighted @ y))
        log_radius[start:start + len(xx)] = -epsilon * log_z
    return normals, log_radius


def refract_using_normals(incident, normals, kappa):
    """Transmit unit incident rays through normals, with kappa=n2/n1.

    Normal magnitude and sign do not matter. Choose the transmitted branch
    on the same normal side as the incident propagation direction. Raise
    on total internal reflection rather than discarding or clipping rays.
    """
    incident = _unit_vectors(incident, "incident")
    normals = np.asarray(normals, dtype=float)
    if normals.shape != incident.shape or not np.isfinite(normals).all():
        raise ValueError("normals must be finite and match incident rays")
    if not np.isfinite(kappa) or kappa <= 0:
        raise ValueError("kappa=n2/n1 must be finite and positive")
    lengths = np.linalg.norm(normals, axis=1)
    if (lengths == 0).any():
        raise ValueError("Surface normals cannot be zero")
    n = normals / lengths[:, None]
    a = np.einsum("ij,ij->i", incident, n)
    n = np.where((a < 0)[:, None], -n, n)
    a = np.clip(np.abs(a), 0, 1)
    discriminant = 1 - (1 - a*a) / (kappa*kappa)
    if (discriminant < -1e-12).any():
        raise ValueError(f"Total internal reflection for {int((discriminant < -1e-12).sum())} rays")
    transmitted = ((incident - a[:, None]*n) / kappa
                   + np.sqrt(np.maximum(discriminant, 0))[:, None]*n)
    return transmitted
