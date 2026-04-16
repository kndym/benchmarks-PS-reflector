"""
refracter: Python implementation of the Point Source Far-Field Refractor/Reflector
Problem using Entropic Optimal Transport (Sinkhorn divergence).

Default use case: refractor (κ = 0.6).  Set κ = 1.0 for the reflector case.
"""

from .qmc import load_main_cloud, load_small_cloud, load_push_cloud
from .distributions import BENCHMARKS, stereo_north, stereo_south
from .cost import cost_matrix_chunk
from .sinkhorn import run_sinkhorn_divergence
from .build import build_reflector, c_transform_gc, c_transform_fc, build_regular_grid, reflector_on_regular_grid
from .pushforward import ray_trace

__all__ = [
    "load_main_cloud",
    "load_small_cloud",
    "load_push_cloud",
    "BENCHMARKS",
    "stereo_north",
    "stereo_south",
    "cost_matrix_chunk",
    "run_sinkhorn_divergence",
    "build_reflector",
    "c_transform_gc",
    "c_transform_fc",
    "build_regular_grid",
    "reflector_on_regular_grid",
    "ray_trace",
]
