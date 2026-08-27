"""Print differences between the chunked cost matrix and a basic loop."""

import math
import os
import sys

# Run without creating or changing Python bytecode files in the repository.
sys.dont_write_bytecode = True

# Make ``python tests/test_cost_matrix.py`` work from the repository root.
repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, repo_root)

import numpy

from refracter.cost import cost_matrix_chunk, set_kappa


EPS_CLIP = 1e-15


def naive_cost_matrix(x, y, kappa):
    """Calculate costs using only Python loops and math.log."""
    matrix = []
    for x_vec in x:
        row = []
        for y_vec in y:
            dot = 0.0
            for i in range(len(x_vec)):
                dot += x_vec[i] * y_vec[i]
            argument = max(1.0 - kappa * dot, EPS_CLIP)
            row.append(-math.log(argument))
        matrix.append(row)
    return matrix


def chunked_cost_matrix(x, y, chunk_size):
    """Calculate the full matrix by calling cost_matrix_chunk in chunks."""
    matrix = []
    for start in range(0, len(x), chunk_size):
        block = cost_matrix_chunk(x[start : start + chunk_size], y)
        for row in block:
            matrix.append([float(value) for value in row])
    return matrix


def random_unit_vectors(count, rng):
    """Make unit vectors; NumPy is used only to generate random numbers."""
    vectors = []
    for _ in range(count):
        vector = [float(rng.normal()), float(rng.normal()), float(rng.normal())]
        length = math.sqrt(sum(value * value for value in vector))
        vectors.append([value / length for value in vector])
    return vectors


rng = numpy.random.default_rng(20260826)
x = random_unit_vectors(1000, rng)
y = random_unit_vectors(1200, rng)


for kappa in (1.0, 0.6):
    set_kappa(kappa)
    expected = naive_cost_matrix(x, y, kappa)

    for chunk_size in (1, 17, 64, 256, len(x) + 1):
        actual = chunked_cost_matrix(x, y, chunk_size)
        max_distance = 0.0
        for i in range(len(x)):
            for j in range(len(y)):
                distance = abs(actual[i][j] - expected[i][j])
                max_distance = max(max_distance, distance)

        print(
            "kappa=", kappa,
            "chunk_size=", chunk_size,
            "max_absolute_distance=", format(max_distance, ".17g"),
        )
