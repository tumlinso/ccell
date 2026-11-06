# Working on preprocess inside ccell

## Purpose and boundaries

The `preprocess` subtree is part of ccell's independent C++20 reference biology
package. Its public functions provide an explicit CPU reference built on
AnnData++. Scanpy supplies pinned numerical comparison behavior; it does not
dictate Python-style call semantics.
Numerical similarity is required; bitwise identity and universal equality of
intermediates are not requirements. Explain consequential differences and
qualify them with scientific quantities and tolerances.

`preprocess` depends on AnnData++. It does not depend on Python, Scanpy at
runtime, Torch, CUDA, or Cellerator. Keep preprocessing, fitting, PCA,
neighbors, GPU execution, and dataset orchestration as separate concerns.

## C++ and numerical style

- Follow STYLE.md and the protected
  [Cellerator C++ coding style](/home/tumlinson/Cellerator/coding_style.md).
  Keep its reference unchanged.
- Use C++20 named modules and .cppm interfaces. Public consumption is
  import preprocess;. Do not require import std; or add a public umbrella
  header.
- Keep exported module interfaces and consumer recompilation small. Put
  implementation-only code in compiled source files.
- Use direct arrays, pointers, extents, strides, sparse offsets and indices.
  Make allocation, conversion, mutation, and materialization explicit.
- Use AnnData++ views directly. Keep shape, layout, lifetime and buffer
  representation with AnnData++; do not maintain a parallel matrix abstraction.
  Preprocess owns numerical operations and their additional contracts.
- Resolve validation and representation dispatch at acquisition boundaries.
  Keep element loops free of type-erased dispatch, repeated validation and
  hidden synchronization. Callers own parallel execution.
- Keep numerical types and accumulation explicit. Preserve the operation's
  declared mathematical method and test against the pinned reference cases.

## Build and validation

Use Linux, Clang 18 with matching clang-scan-deps, CMake 3.28 or newer, and
Ninja 1.11 or newer.

The maintained source is `ccell/preprocess`. Build and validate it through the
root ccell CMake project as described in [the root README](../README.md). It
has no standalone configure or package workflow. `CCELL_BUILD_TESTS` includes
the prefixed native suites and example; `CCELL_BUILD_QUALIFICATION` enables the
optional Scanpy comparison and native benchmark sources.

Optional Python interoperability is development-only. Its reference environment
pins Scanpy 1.11.5 and the packages in tests/requirements-reference.txt.
See docs/development.md and docs/validation.md.

Performance probes should be bounded to 30 seconds and compare identical
workloads. Keep source and toolchain identity with results. Do not claim speed
or memory improvements from source inspection alone.

## Ownership and preservation

Preserve unknown work and supplied inputs. Do not reset or clean Git state,
force-push, replace authority data, or edit the sibling AnnData++ or Velocity++
components as part of preprocessing work. Keep generated builds, local
reference environments, data and raw reports out of Git.

The former `preprocess-cppm` GitHub repository and Project Control registration
describe historical standalone development. They do not establish current
ownership or workflow authority; the maintained source is `ccell/preprocess`.

Use the configured Codex agent profiles for useful independent work. Root
retains user intent, interface decisions, integration, and acceptance. Use
Codex subagents for engineering; preserve concurrent edits and never use
Project Control local coding workers as engineering delegates.
