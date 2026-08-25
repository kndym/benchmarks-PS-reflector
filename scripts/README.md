# Scripts

The two result-generation commands are thin entry points around
`refracter.refraction_pipeline`:

- `generate_results.py` runs the standard uniform-to-uniform refraction case.
- `generate_results_all_pairs.py` runs all 16 source/target density pairs.

`visualize.py` renders plots from a reflector output directory produced by the
top-level `benchmark_reflector.py` command.

The `reflector/` subdirectory contains optional C++-parity experiments for
the old `output_cpp_*` and `output_py_*` formats. They are not required by
the primary refraction pipeline.
