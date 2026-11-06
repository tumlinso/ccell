# Coding style

Write compact and explicit numerical code in the tradition of classic C and
Fortran scientific libraries. Keep the calculation, storage layout, and
execution visible.

- Keep related mathematics and invariants together. Use concise conventional
  notation and direct loops; extract helpers for real reusable calculations
  or distinct contracts.
- Use arrays, pointers, references, dimensions, and strides that show where
  values come from and where they go. Describe ownership, aliasing, and byte
  versus element strides precisely.
- The standard reference allocates nothing in sizing, validation, preparation,
  binding, or execution. All outputs and scratch are supplied by the caller.
- Use small structs and ordinary free functions for data and computations.
  Introduce classes for ownership or real invariants. Do not add callback
  layers, wrapper hierarchies, or generic frameworks without a concrete user.
- Validate types, layouts, ownership, and dimensions at operation boundaries.
  Keep acquired numeric loops free of repeated checks, type erasure and
  library synchronization.
- Keep C++20 named-module interfaces small and independently compiled.
  Benchmark runtime or build-cost claims on identical workloads.

Follow the binding [ccell numerical style](../../coding_style.md), including its
minimal include policy, strict numerical rules, and downstream replacement
contracts. The Cellerator style is its baseline and remains unchanged. Fast
imports and consumer recompilation are required.
