# Development

Preprocessing is maintained at `src/preprocess/`, with native tests and the
interoperability driver under `tests/preprocess/`. The ccell root build owns the
module target, tests, installation, and optional qualification targets. Use the
Clang, clang-scan-deps, CMake, and Ninja versions in the
[root build instructions](../../README.md).

The numerical module exposes compact operation contracts over acquired
AnnData++ views and caller supplied outputs and scratch. Its standard path
allocates no memory during sizing, validation, preparation, or execution.
Consumers create AnnData++ storage, retain its views, allocate and initialize
all output and workspace, and control parallelism. Checked calls report status
and required capacity; trusted calls use the documented preconditions.

With the reference Python environment installed, run the independent Scanpy
1.11.5 comparison from the repository root:

~~~sh
python -m pip install -r tests/preprocess/requirements-reference.txt
cmake --preset debug -DCCELL_BUILD_QUALIFICATION=ON
cmake --build build/debug --target ccell_preprocess_interop_driver --parallel 2
python tools/preprocess/qualify.py --driver build/debug/preprocess/ccell_preprocess_interop_driver
~~~

By default, the report is written to
`build/qualification/preprocess-reference.json`; use `--output` for another
disposable path. Historical reports under `tests/preprocess/results/` are
inputs to provenance and must not be overwritten.

The native runtime observation uses the same root option:

~~~sh
cmake --build build/debug --target ccell_preprocess_benchmark --parallel 2
python tools/preprocess/measure_runtime.py --driver build/debug/preprocess/ccell_preprocess_benchmark
~~~

Its report defaults to `build/qualification/preprocess-runtime.json`. Each
subprocess is bounded to 30 seconds; reports characterize the stated inputs
and toolchain and do not establish universal speed or biological validity.

Older standalone configure, package, qualification, and build measurements
remain historical evidence. Root CMake and the current source tree are the
authority for maintained builds.
