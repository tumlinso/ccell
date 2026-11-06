# preprocess

This source is maintained inside ccell at `preprocess/` as part of ccell's
independent C++20 reference biology package. It is ordinary ccell source, not a
separately maintained GitHub project or standalone CMake package.

The `preprocess` C++20 module provides reusable preprocessing calculations for
annotated biological data. It requires
[AnnData++](https://github.com/tumlinso/anndatapp) and provides a direct
CPU reference based on Scanpy's numerical behavior.

The interface is C++-native: typed options and results, explicit ownership,
mutable versus read-only array access, and free functions. Cellerator may use
the calculations or compare specialized implementations with the reference.
The package has no runtime dependency on Python, Scanpy, Torch, CUDA, or
Cellerator.

The public named module remains `preprocess`; ccell also re-exports it through
`import ccell.preprocess;`.

The example below assumes an existing `anndatapp::AnnData` named `data`; a
complete compiled example is in [examples/basic.cpp](examples/basic.cpp).

~~~cpp
#include <utility>
import preprocess;

auto result = preprocess::qc(preprocess::source(data));
auto cells = preprocess::filter_cells(result, {.min_detected = 200});

auto normalized = preprocess::floating_copy(
    preprocess::source(data), anndatapp::DType::float32);
auto factors = preprocess::normalize_total(
    normalized.edit(), {.target_sum = 1e4});
preprocess::log1p(normalized.edit());
preprocess::annotate_qc(data, std::move(result));
preprocess::annotate_normalization(data, std::move(factors));
data.set_layer("normalized", std::move(normalized));
~~~

The numerical scope includes quality-control summaries, filter masks, total
normalization, log1p, feature scaling, and Seurat or Cell Ranger dispersion
based highly variable gene selection. It does not include PCA, neighborhood
construction, Seurat v3, dataset orchestration, or device execution. See
[operations](docs/operations.md) for numerical conventions and
[development](docs/development.md) for builds and qualification.

## Build and validation

Build, test, and install this component through the ccell root project. Follow
the [root build instructions](../README.md); there is no component-level
configure step or `find_package(preprocess)` package. The root build exposes
the module target `ccell::preprocess` and the component module import above.
Set `CCELL_BUILD_TESTS=ON` to include the native preprocessing suites and
example. Set `CCELL_BUILD_QUALIFICATION=ON` to build the optional Scanpy
comparison driver and native benchmark.

The former `preprocess-cppm` GitHub repository and Project Control registration
are historical records, not current source ownership or workflow authority.

The prior standalone Scanpy qualification recorded 302 cases across supported
operations, layouts, precisions, annotations, and declared numerical
conventions. See [validation](docs/validation.md) for that historical evidence
and its coverage limits.

## License

Original code is MIT licensed; see the ccell [LICENSE](../LICENSE) and
[third-party notices](../THIRD_PARTY_NOTICES.md).
