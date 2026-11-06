# Allocation-free reference cutover

The reference uses caller-supplied inputs, outputs, preparation storage and
scratch. AnnData++ acquisition and ownership belong to the caller. Downstream
replacements may allocate explicitly and use different mathematics and layouts.
The binding source rules are in [coding_style.md](../coding_style.md).

## Current checks

On Linux, Clang 18 with matching clang-scan-deps, CMake 3.28.3 and Ninja:

- Root Debug, Release, shared, and ASan/UBSan configurations each passed 10/10
  CTests, including umbrella/component imports and installed consumers.
- The root sanitizer build instruments AnnData++, numerical libraries, and
  consumers. Leak detection was disabled (`ASAN_OPTIONS=detect_leaks=0`).
- Standalone Velocity++ built all targets, including seven qualification
  drivers/benchmarks, and passed 14/14 CTests including its installed consumer.
- The shipped-source audit passed 29 C++ files. AnnData++ is excluded as the
  explicit data-model exception. Development-only harnesses may use STL.
- Ordinary Debug/Release probes observed zero public C allocator and C++ new
  calls across covered reference phases, fresh/repeated and error/capacity
  paths. Private allocators and direct memory mapping are outside interception.
  Under ASan the incompatible allocator interposition is disabled explicitly;
  the same numerical cases run as `ccell-preprocess-sanitized-behavior`.
- Targeted independent review findings were repaired and re-reviewed: capacities,
  empty reductions, sparse traversal, writable overlap/orientation, scientific
  phase composition, signed strides, and partial-output invalidation.

Logs and source-bound reports are under `build/refactor/`. No runtime speedup
or biological-validity claim follows from these checks.

## Independent numerical comparisons

Fresh baseline snapshots preceded numerical edits. Current reports and matching
before/after source/binary manifests are in `build/refactor/current/`:

- Scanpy 1.11.5: 302/302 preprocessing cases passed with frozen source/driver
  hashes. Existing tolerances, feature alignment, ties, layouts and precision
  conventions were retained.
- Pinned scVelo deterministic comparison: 4/4 passed.
- Stochastic estimator, independent FP64 initial-fit, and production-class
  gates: 32/32 passed. The exact same four baseline initial-offset source
  diagnostics remain; native matches the independent FP64 calculation. These
  source misses are recorded as misses, with no widened tolerances.
- Sensitivity: the existing single-fixture check passed; reported quantities
  match baseline apart from the rebuilt driver hash.

Historical receipts remain unchanged. Qualifiers and measurement tools now
write new results to build-owned paths rather than historical receipt files.

## Module rebuild contract

`tools/check_module_rebuild.py` built a disposable full-source copy in Debug
and Release, each passing 10/10 tests. A no-op preserved 106 compile artifacts.
Consumer edits rebuilt only affected root/installed consumer objects. Private
preprocessing and Velocity++ edits rebuilt their implementation objects while
preserving module files, module-interface objects and consumer objects by hash
and timestamp. A separate exported-specialization edit rebuilt its module and
importer and changed the executable result as expected.

The source manifest and complete stamps are in `build/refactor/module-rebuild.json`
and `build/artifacts/module-rebuild/20261106T150327.221480Z/`. Later changes were
confined to development tools/tests/docs and sanitizer test registration;
module interfaces and numerical implementation sources matched the snapshot.

AnnData++ remains at e8c51ce7c3d8c2ded9d1b5f43c1e4d035ea13629 with no source
edits. Its Git identity, Velocity++ history, preprocessing recovery records,
and historical scientific fixtures/receipts are preserved.
