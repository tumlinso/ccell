# ccell

A C++20 biology package with explicit reference calculations, native data
access, and composable named modules. AnnData++ provides the data model,
preprocessing is local source, and Velocity++ provides RNA velocity.
Cellerator can consume these reference contracts and supply specialized
implementations.

```cpp
import ccell;             // Re-exports all three components.

// Or import components individually:
import ccell.anndata;     // Re-exports anndatapp.
import ccell.preprocess;  // Re-exports the local preprocess module.
import ccell.velocity;    // Re-exports velocitypp and velocitypp.precision.
```

The component imports expose `ccell::anndata`, `ccell::preprocess`, and
`ccell::velocity` as namespace aliases. The original namespaces and standalone
AnnData++ and Velocity++ modules remain available. These imports share the
same types, ownership, and implementations; they do not copy or wrap datasets.
The dots are part of C++ named-module names, not a separate submodule mechanism.

## Build

Use Clang 18 or newer with matching `clang-scan-deps`, CMake 3.28 or newer,
and Ninja 1.11 or newer. Initialize the AnnData++ and Velocity++ submodules:

```sh
git submodule update --init --recursive
cmake --preset debug
cmake --build --preset debug
ctest --preset debug
```

The debug and release presets build without storage backends. H5AD and Zarr
can be enabled through `ANNDATAPP_ENABLE_HDF5` and `ANNDATAPP_ENABLE_ZARR`;
see [AnnData++ storage dependencies](anndata/docs/storage-dependencies.md).
AnnData++'s own defaults apply when configuring without these presets.
Python, Torch, CUDA, and Cellerator are not runtime or build dependencies of
the reference calculation package. Velocity++ accepts FP32/FP64 arrays;
the FP16 implementation remains in the separate Cellerator checkout.

The root build compiles the component sources together. It does not require
preinstalled copies of AnnData++, Velocity++, or preprocessing.

## Installed consumers

```sh
cmake --install build/debug --prefix /path/to/ccell-install
```

```cmake
cmake_minimum_required(VERSION 3.28)
project(example LANGUAGES CXX)
set(CMAKE_CXX_EXTENSIONS OFF)
find_package(ccell CONFIG REQUIRED)
add_executable(example main.cpp)
target_link_libraries(example PRIVATE ccell::ccell)
```

Configure the consumer with the supported compiler, Ninja, and
`-DCMAKE_PREFIX_PATH=/path/to/ccell-install`. The installation includes the
component libraries, module sources, and CMake metadata. Consumers compile
module interfaces with their own compatible toolchain.
Preprocessing is part of the `ccell` build and installation; it has no separate
configure step or CMake package.

## Source layout

- `anndata/`: AnnData++ submodule; annotated data and optional storage adapters.
- `velocity/`: Velocity++ submodule; standalone FP32/FP64 reference algorithms.
- `preprocess/`: ordinary ccell source; QC, filtering, normalization, scaling,
  and highly variable gene selection.
- `src/ccell.*.cppm` and `ccell.ccm`: component re-exports and umbrella module.

The Velocity++ submodule uses
[tumlinso/velocitypp](https://github.com/tumlinso/velocitypp), with the local
source directory retained as `velocity/`. Its C++ module, namespace, and CMake
package use `velocitypp`.

The preprocessing repository history is preserved locally in
`.git/modules/preprocess`; `.ccell-local/` records its original Git pointer
and revision. The preprocessing files are ordinary ccell source rather than
a submodule. No remote repository was deleted.

The root [LICENSE](LICENSE) retains the preprocessing implementation's MIT
notice. [Third-party notices](THIRD_PARTY_NOTICES.md) retain Scanpy's notice and
the bundled dependencies' notices. Both files are included in installations.
