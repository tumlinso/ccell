# ccell

ccell is a C++20 numerical development reference for single-cell methods. It
exposes calculations, arrays, pointers, layouts, working storage, and execution
to competent developers. The standard reference allocates nothing: callers
provide inputs, outputs, preparation storage, intermediate buffers, and scratch.
AnnData++ owns biological data and storage; callers acquire it explicitly before
reference calculation. Caller threads own scheduling and synchronization.

Consumers may replace a primitive, connected region, or complete algorithm with
novel mathematics and representations. Compile-time templates and runtime typed
function pointers support caller-defined state. Custom layouts can persist
through an entire replacement region; conversion is required only where that
region connects to retained code. Replacements may allocate explicitly.
Cellerator is a downstream consumer, not a dependency of this package.

Read the binding [coding style](coding_style.md) before numerical development.
It extends Cellerator's compact numerical style with minimal scalar includes,
compiler builtins, visible data structures, and allocation-free reference work.
Development-only tests and qualification tools may use ordinary library code.
AnnData++ is the explicit exception for the existing data model and ownership.

```cpp
import ccell;             // All three components, with the same component types.
// Or:
import ccell.anndata;     // Re-exports anndatapp.
import ccell.preprocess;  // Re-exports preprocess.
import ccell.velocity;    // Re-exports velocitypp and velocitypp.precision.
```

The namespace aliases are `ccell::anndata`, `ccell::preprocess`, and
`ccell::velocity`. Original module and namespace identities remain available.
The component imports do not wrap or copy datasets. Existing STL-bearing and
allocating numerical calling forms are replaced outright; see the component
API documentation for caller-storage interfaces.

## Build and install

Use Clang 18 with matching `clang-scan-deps`, CMake 3.28 or newer, and Ninja:

```sh
git submodule update --init --recursive
cmake --preset debug
cmake --build --preset debug
ctest --preset debug
cmake --install build/debug --prefix /path/to/ccell-install
```

Release, shared, and sanitizer presets are also provided. Storage backends are
off in these presets; enable `ANNDATAPP_ENABLE_HDF5` or `ANNDATAPP_ENABLE_ZARR`
explicitly when wanted. [Storage dependencies](external/anndatapp/docs/storage-dependencies.md)
remain private to AnnData++. Python, Scanpy, Torch, CUDA, and Cellerator are not
runtime dependencies of the reference. Numerical storage is FP32/FP64 where
supported; Cellerator owns specialized device and FP16 implementations.

Installed CMake consumers use `find_package(ccell CONFIG REQUIRED)` and link
`ccell::ccell`. Configure with the compatible compiler, Ninja, and
`-DCMAKE_PREFIX_PATH=/path/to/ccell-install`. Module sources and metadata are
installed so consumers compile interfaces with their own compatible toolchain.
Ordinary calculations are separately compiled; selective template customization
has its explicit consumer compilation cost. Preprocessing has no standalone
build or installed CONFIG package.

## Ownership and source layout

- `src/ccell.cppm` and component interfaces: umbrella and type-preserving exports.
- `src/preprocess/`: ccell-owned preprocessing calculations.
- `external/anndatapp/`: unchanged AnnData++ Git submodule.
- `external/velocitypp/`: standalone Velocity++ Git submodule.
- Root tests, examples, tools, docs, benchmarks, and planning records group
  preprocessing material under their component directories.

Git submodule identities and independent histories are preserved. The former
preprocessing repository history remains locally in `.git/modules/preprocess`;
`.ccell-local/` retains its recovery records. Historical scientific receipts
remain historical and do not qualify a changed source snapshot.

The root [LICENSE](LICENSE) and [third-party notices](THIRD_PARTY_NOTICES.md)
retain original and upstream notices and are installed with the package.
