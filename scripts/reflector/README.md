# Reflector comparison helpers

These scripts compare the Python reflector experiments with C++ reference
outputs. The supported reflector benchmark entry point is the top-level
`benchmark_reflector.py`; these helpers are only needed when reproducing the
legacy `output_cpp_*`, `output_py_*`, or `output_python/` comparisons.

The scripts expect to be run from the repository root, for example:

```bash
python scripts/reflector/run_compare.py 1600
python scripts/reflector/compare_nk.py 1600
```
