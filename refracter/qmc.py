"""
refracter/qmc.py

Parse and cache QMC point clouds from C++ header files.
"""

import re
import os
import numpy as np

# ---------------------------------------------------------------------------
# Paths to the header files (relative to this file's location)
# ---------------------------------------------------------------------------
_HERE = os.path.dirname(os.path.abspath(__file__))
_BENCHMARK_ROOT = os.path.join(_HERE, "..", "BenchmarkCode")

_MAIN_H   = os.path.join(_BENCHMARK_ROOT, "QuasiMonteCarlo", "MonteCarlo_Pointcloud_3D_128.h")
_SMALL_H  = os.path.join(_BENCHMARK_ROOT, "SmallGrid", "3D_MonteCarlo_Pointcloud_small.h")
_PUSH_H   = os.path.join(_BENCHMARK_ROOT, "PushForward", "PushForward_Cloud_128.h")

# Expected sizes
_NK       = 16488
_NK_SMALL = 381
_DIM      = 3

# ---------------------------------------------------------------------------
# Module-level caches
# ---------------------------------------------------------------------------
_main_cloud_cache  = None
_small_cloud_cache = None
_push_cloud_cache  = None


def _parse_h_array(text: str, varname: str, n_rows: int, n_cols: int = 3) -> np.ndarray:
    """Extract a 2-D numeric array from a C++ header file.

    The function looks for the declaration

        <type> <varname>[...][...] = { ... };

    and then collects all floating-point numbers inside the braces using
    a regex that matches C++ scientific notation (e.g. ``-6.787942e-01``).

    Parameters
    ----------
    text : str
        Full text of the header file.
    varname : str
        C variable name to search for.
    n_rows : int
        Expected number of rows in the array.
    n_cols : int, optional
        Number of columns (default 3).

    Returns
    -------
    np.ndarray of shape (n_rows, n_cols), dtype float64.
    """
    # Find the opening brace that belongs to this variable declaration.
    # We search for "<varname>..." followed by "= {" possibly with whitespace.
    pattern = re.compile(
        r'\b' + re.escape(varname) + r'\s*\[.*?\]\s*(?:\[.*?\]\s*)?=\s*\{',
        re.DOTALL
    )
    m = pattern.search(text)
    if m is None:
        raise ValueError(
            f"Could not find variable '{varname}' in the provided header text."
        )

    start = m.end()  # position right after the opening '{'

    # Find the matching closing brace.  We need to track brace depth because
    # the initializer list might contain nested braces for row delimiters.
    depth = 1
    pos = start
    while pos < len(text) and depth > 0:
        ch = text[pos]
        if ch == '{':
            depth += 1
        elif ch == '}':
            depth -= 1
        pos += 1

    body = text[start: pos - 1]  # exclude the final '}'

    # Extract all numbers.  The pattern matches optional sign, digits, optional
    # decimal point, more digits, and an exponent part — exactly how C++ emits
    # floating-point literals in scientific notation.
    numbers = re.findall(r'[+-]?\d*\.?\d+[eE][+-]?\d+', body)

    total = n_rows * n_cols
    if len(numbers) < total:
        raise ValueError(
            f"Expected at least {total} numbers for '{varname}' "
            f"({n_rows}×{n_cols}), but found only {len(numbers)}."
        )

    arr = np.array([float(v) for v in numbers[:total]], dtype=np.float64)
    return arr.reshape(n_rows, n_cols)


def load_main_cloud():
    """Load the main QMC point cloud (NK=16488 points, dim=3).

    Returns
    -------
    x : np.ndarray, shape (16488, 3)
        Source points on the upper hemisphere.
    y : np.ndarray, shape (16488, 3)
        Target points on the lower hemisphere.
    """
    global _main_cloud_cache
    if _main_cloud_cache is not None:
        return _main_cloud_cache

    with open(_MAIN_H, "r", encoding="utf-8") as fh:
        text = fh.read()

    x = _parse_h_array(text, "x", _NK, _DIM)
    y = _parse_h_array(text, "y", _NK, _DIM)

    _main_cloud_cache = (x, y)
    return _main_cloud_cache


def load_small_cloud():
    """Load the small warm-start QMC point cloud (NK_small=381 points).

    Returns
    -------
    x_small : np.ndarray, shape (381, 3)
    y_small : np.ndarray, shape (381, 3)
    """
    global _small_cloud_cache
    if _small_cloud_cache is not None:
        return _small_cloud_cache

    with open(_SMALL_H, "r", encoding="utf-8") as fh:
        text = fh.read()

    x_small = _parse_h_array(text, "x_small", _NK_SMALL, _DIM)
    y_small = _parse_h_array(text, "y_small", _NK_SMALL, _DIM)

    # The header comment states: "where reflector cost is applied, last values
    # of y should change sign."  The y_small array is stored with the same
    # x-coordinates as x_small but needs z negated to place the points on the
    # lower hemisphere (matching the sign convention of the main y cloud).
    y_small = y_small.copy()
    y_small[:, 2] = -y_small[:, 2]

    _small_cloud_cache = (x_small, y_small)
    return _small_cloud_cache


def load_push_cloud():
    """Load the push-forward point cloud (16488 points on the upper hemisphere).

    Returns
    -------
    push_cloud : np.ndarray, shape (16488, 3)
    """
    global _push_cloud_cache
    if _push_cloud_cache is not None:
        return _push_cloud_cache

    with open(_PUSH_H, "r", encoding="utf-8") as fh:
        text = fh.read()

    # The variable is named Push_Cloud and has size Push_Cloud_Size (=16488)
    push_cloud = _parse_h_array(text, "Push_Cloud", _NK, _DIM)

    _push_cloud_cache = push_cloud
    return _push_cloud_cache
