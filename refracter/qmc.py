# Parse and cache quasi-Monte Carlo point clouds from C++ headers.

import re
import os
import numpy as np

# Resolve the headers relative to this module, not the current directory.
_HERE = os.path.dirname(os.path.abspath(__file__))
_BENCHMARK_ROOT = os.path.join(_HERE, "..", "BenchmarkCode")

_MAIN_H   = os.path.join(_BENCHMARK_ROOT, "QuasiMonteCarlo", "MonteCarlo_Pointcloud_3D_128.h")
_SMALL_H  = os.path.join(_BENCHMARK_ROOT, "SmallGrid", "3D_MonteCarlo_Pointcloud_small.h")
_PUSH_H   = os.path.join(_BENCHMARK_ROOT, "PushForward", "PushForward_Cloud_128.h")

# Expected sizes in the reference point-cloud files.
_NK       = 16488
_NK_SMALL = 381
_DIM      = 3

# Cache parsed arrays so repeated benchmark stages do not reread the headers.
_main_cloud_cache  = None
_small_cloud_cache = None
_push_cloud_cache  = None


def _parse_h_array(text: str, varname: str, n_rows: int, n_cols: int = 3) -> np.ndarray:
    # Find a named C++ initializer and return its first n_rows*n_cols numbers.
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
    # Load the 16,488-point source and target clouds used by the reflector case.
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
    # Load the 381-point clouds used for the reflector warm start.
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
    # Load the source cloud used to test the reflector push-forward.
    global _push_cloud_cache
    if _push_cloud_cache is not None:
        return _push_cloud_cache

    with open(_PUSH_H, "r", encoding="utf-8") as fh:
        text = fh.read()

    # The variable is named Push_Cloud and has size Push_Cloud_Size (=16488)
    push_cloud = _parse_h_array(text, "Push_Cloud", _NK, _DIM)

    _push_cloud_cache = push_cloud
    return _push_cloud_cache
