# ccell preprocessing

Preprocessing is ccell owned source under `src/preprocess/`. It provides direct
C++20 numerical operations over AnnData++ arrays for single cell development.
The standard implementation makes traversal explicit and performs no
allocation. Callers acquire storage, supply outputs and scratch, and choose
threading and synchronization. See the [numerical contract](operations.md).

The `preprocess` named module is available directly and through
`import ccell.preprocess;`. AnnData++ remains the owner of array storage and
views. Preprocessing does not introduce a parallel matrix or ownership model.
The API exposes reference operations and permits replacement of individual
operations or larger regions with caller defined algorithms and layouts.

The scope includes quality control and filters, feature statistics, total
normalization, log transforms, scaling, and Seurat or Cell Ranger highly
variable gene calculations. The package has no runtime dependency on Python,
Scanpy, Torch, CUDA, or Cellerator. Scanpy 1.11.5 is used only for optional
independent numerical qualification.

## Build and validate

Build and install preprocessing through the ccell root CMake project. Follow
the [root build instructions](../../README.md). The component has no standalone
package or configure step. The root build provides the `ccell::preprocess`
target and optional native tests, example, Scanpy comparison driver, and
runtime benchmark.

The historical Scanpy report records 302 cases for an earlier source snapshot.
It is preserved as history and does not qualify the relocated or current code.
See [development](development.md) for current commands and
[validation](validation.md) for retained fixtures and evidence.

The original code is MIT licensed. See the ccell root [LICENSE](../../LICENSE)
and [third-party notices](../../THIRD_PARTY_NOTICES.md).
