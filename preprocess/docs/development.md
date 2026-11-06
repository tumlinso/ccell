# Development

## Current build

Preprocessing is ordinary source within ccell and builds through the root CMake
project. Follow the [root build instructions](../../README.md) for toolchain,
configuration, testing, and installation. The root build provides the stable
target `ccell::preprocess` and re-exports this component through
`import ccell.preprocess;`; the original `import preprocess;` module remains
available. There is no standalone preprocess project or installed CMake package.

With `CCELL_BUILD_TESTS=ON`, the root build includes
`ccell_preprocess_numerics`, `ccell_preprocess_hvg`, and
`ccell_preprocess_example`. The optional Scanpy comparison driver and native
benchmark are enabled by the root option `CCELL_BUILD_QUALIFICATION`.

## Optional numerical qualification

Install the pinned development dependencies from
`preprocess/tests/requirements-reference.txt` in an isolated environment. From
the ccell root, build and run the independent Scanpy 1.11.5 comparison:

~~~sh
python -m pip install -r preprocess/tests/requirements-reference.txt
cmake --preset debug -DCCELL_BUILD_QUALIFICATION=ON
cmake --build build/debug --target ccell_preprocess_interop_driver --parallel 2
python preprocess/tools/qualify.py --driver build/debug/preprocess/ccell_preprocess_interop_driver
~~~

The native runtime probe uses the same root option and benchmark target. Choose
an allowed idle CPU and keep the output under this component's results folder:

~~~sh
cmake --preset debug -DCCELL_BUILD_QUALIFICATION=ON
cmake --build build/debug --target ccell_preprocess_benchmark --parallel 2
python3 preprocess/tools/measure_runtime.py --driver build/debug/preprocess/ccell_preprocess_benchmark --cpu 2 --output preprocess/benchmarks/results/runtime.json
~~~

Runtime observations are medians of three fresh processes. Each measured
subprocess is capped at 30 seconds. They characterize the supplied workload;
they are not comparisons or speedup claims.

## Historical standalone workflow

The original v0.1 qualification ran when preprocess was a standalone CMake
package. Its helper and presets have been retired; these commands are retained
only as history and are not current build instructions:

~~~sh
sh tools/build_anndatapp.sh /path/to/anndatapp
cmake --workflow --preset native-check
~~~

The corresponding source, toolchain, numerical, and measurement evidence is
preserved in [validation](validation.md) and the result files. The former
`preprocess-cppm` GitHub repository and Project Control registration record
that standalone period; they do not establish current ownership or workflow
authority.
