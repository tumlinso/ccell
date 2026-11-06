# ccell numerical development style

ccell is a C++20 numerical development reference for competent programmers who
require direct, potentially unsafe control of biological data, memory,
algorithms, representations, and execution. Expose the machine and calculation.
The reference is usable code and scientific comparison material; its internal
architecture is not a requirement imposed on downstream implementations.

## Authority and scope

Read this note before numerical design, implementation, or review. It is binding
contributor authority. Do not weaken, replace, move, or remove it, or weaken its
protection in AGENTS.md, without explicit user authorization for that change.
General refactoring authority does not authorize such changes.

The compact mathematical style of Cellerator's protected coding_style.md is the
baseline, strengthened here. Preserve that independent document unchanged.
This note is self-contained: a checkout does not need Cellerator to build.

These rules cover shipped implementations, public templates, and examples,
including Velocity++. AnnData++ is the explicit data-model/ownership exception;
its imported implementation is not permission to introduce library machinery
into numerical code. Development-only tests, qualification drivers, Python
scripts, and measurement infrastructure may use ordinary library facilities.
Their dependencies must not leak into shipped interfaces or implementations.

## Write the calculation directly

Use compact C-style numerical code: pointers, arrays, indices, extents, strides,
loops, arithmetic, plain structs, and deliberate function calls. Keep related
mathematics and its invariants together. Source must reveal traversal and data
layout independently of whether a compiler makes it fast.

Use short names tied to established notation: i, j, k, m, n, lda, dt, p. Give
public operations enough semantic context to stand alone. Show consequential
storage, multiplication, accumulation, and output types and their conversions.
Use auto deliberately to distinguish values, aliases, and pointers.

Keep mathematical expressions whole, such as y[i] = a*x[i] + b*y[i]. Introduce
intermediates for reuse, mathematical significance, or evaluation semantics.
Use ternaries for value selection and compact piecewise definitions; use direct
control flow for distinct actions. Preserve types, rounding, sequencing,
aliasing, and lifetimes when condensing an expression.

Prefer compact guards, loops, and small functions. Group related declarations
and statements. Keep coherent algorithms together; there is no function-length
quota. Extract a function for a real computation, reusable operation, or distinct
contract, rather than to meet a line count. Deliberate duplication is reasonable
when numerical treatment, traversal, or ownership differs.

## Data and language facilities

Use C++ language features to expose useful structure: templates, references,
enums, constexpr, if constexpr, small callable structs, and local lambdas.
Compile-time machinery must remove identifiable repetition or recurring work.
Avoid generic frameworks, registration systems, class hierarchies, forwarding
layers, and policy machinery without a concrete computational need.

No STL containers or algorithms in shipped numerical code: no vector, list,
map, set, string, optional, span, ranges, function, smart pointers, or equivalent
machinery concealed behind another name. Do not replace vector with a home-made
growing container. Operation-specific state consists of visible fields, pointers,
counts, capacities, and explicit initialization. Classes require real invariants.

AnnData++ owns the existing array/view/lifetime model. Reuse its acquired
descriptors where they express the contract; do not create a second general
matrix or storage model. Operation-specific topology and workspace are ordinary
borrowed algorithm state, not another ownership framework.

## Includes and compiled modules

Every include requires the declarations it supplies. Use the smallest necessary
fundamental headers: stddef.h, stdint.h, and float.h. Prefer compiler math and
memory builtins over broad runtime headers. Preserve their required runtime
linkage and numerical semantics; this is not a freestanding-runtime promise.
Do not include cstdlib, stdlib.h, cmath, math.h, algorithm, numeric, utility,
type_traits, limits, or other broad C++ library headers in numerical code.
Audit transitive includes originating in authored code as well as direct ones.

Keep public named-module interfaces small. Ordinary calculations and common
specializations belong in separately compiled implementation units. Export
necessary declarations/descriptors and small, genuine composition templates.
Make wider compile-time customization an explicit consumer choice. An
implementation-only edit must leave public module artifacts and consumer
objects unchanged; relinking the changed library is expected.

## The standard reference allocates nothing

Sizing, validation, preparation, binding, and execution use caller-supplied
storage. No new/delete, malloc/free, growing buffers, owning result constructors,
materialization, hidden caches, or allocation-capable acquisition inside the
standard reference. This includes auxiliary diagnostics, topology, selection,
sensitivity, and error reporting. Empty cases obey the same rule.

Report sufficient capacities and alignment before work. Accept explicit input,
output, intermediate, and scratch storage, report used lengths, and reuse it.
Distinguish byte counts from element counts. Initialization and zeroing are
visible operations. Capacity failure reports required capacity rather than
allocating or growing storage. No hidden rollback or full-data staging.

The caller performs AnnData++ storage creation, metadata attachment, selection
materialization, and other potentially allocating acquisition. An already
acquired retained descriptor may be copied only when that path allocates
nothing. Fully borrowed pointer calls leave lifetime entirely with the caller.

## Contracts, unsafe access, and execution

Make mutability, pointer stability, lifetime, shape, signed strides, sparse
index/offset widths, overlap, alignment, and invalidation explicit. Expose
unchecked entry points. Validate once through explicit operations and let trusted
numerical routines rely on their preconditions. Do not repeat acquisition,
type-erased dispatch, validation scans, or synchronization in element loops.

Reference validation reports compact status, location, and required capacity;
it does not construct allocating exception messages. Report scientific
outcomes such as singular fits separately from argument/capacity failures.
Document partial-output validity and completed work when execution can fail.
Unsafe access transfers stated obligations; it does not change C++ semantics.

Callers own threads, scheduling, worker-local storage, reduction order, and
synchronization. Expose useful blocks, partial reductions, and merge operations.
Reference CPU calls complete synchronously. Do not introduce background pools,
mutable global dispatch registries, implicit synchronization, or hidden state.

Ordinary labels and goto are permitted when they expose phases, retry/fallback,
cleanup, or exits from nested loops. Computed goto is permitted as a documented
compiler extension inside compiled implementation functions. Keep targets in
the same function and obey lifetime/destruction restrictions. Label addresses
are not an external callback or persistent identity contract.

Alignment/nonalias assumptions and target-specific intrinsics are permitted
with precise preconditions. Preserve declared numerical precision, evaluation,
convergence, ordering, and fallback contracts. No blanket fast-math changes.

For example, a direct numerical entry point exposes both its writes and capacity:

```cpp
int row_sums(const double* x, size_t m, size_t n, ptrdiff_t sr, ptrdiff_t sc,
             double* y, size_t capacity) {
    if (capacity < m) return 1;
    for (size_t i = 0; i < m; ++i) {
        double s = 0;
        for (size_t j = 0; j < n; ++j) s += x[ptrdiff_t(i)*sr + ptrdiff_t(j)*sc];
        y[i] = s;
    }
    return 0;
}
```

The caller supplies valid storage, non-overlap, and signed element strides.

## Downstream replacement freedom

Consumers may replace primitives, connected regions, fused blocks, or entire
algorithms with novel methods, layouts, intermediates, and execution strategies.
Only connections to retained code impose compatibility obligations. No consumer
must reproduce reference phases, host arrays, workspace, or intermediate types.
Custom representations may persist across the entire replacement region.

Expose meaningful free subroutines and narrow composition points. Templates
support arbitrary concrete consumer state and results. Runtime bindings use
typed function pointers and caller-owned context; custom data and results stay
in that context. Consumers synchronize binding changes between invocations and
own any device completion/dependency contract. Avoid per-element dispatch.

Replacements may allocate where necessary. Allocation amount, reason, lifetime,
copying, and synchronization must remain visible; avoid incidental allocation.
This exception does not permit allocation inside standard reference routines.
Style compliance and scientific equivalence are separate questions. A novel
method declares its own mathematical contract; an equivalence claim needs
appropriate evidence rather than resemblance to reference control flow.

A runtime replacement can keep its entire representation in caller state:

```cpp
struct Phase { void* state; int (*run)(void*); };
int execute(Phase p) { return p.run(p.state); }

template<class Operation> auto execute(Operation& op) { return op.run(); }
```

The concrete callback interprets its own state. Neither binding requires a
reference matrix, intermediate, or result type. The runtime phase is bound
before execution; a consumer may choose a completely different function/state
pair on its next call. These small bindings are examples, not a required generic
execution framework.

## Comments, review, and evidence

Whitespace should reveal terms, factors, phases, and loop structure. Comments
supply units, ordering, ownership, lifetimes, preconditions, numerical meaning,
and measured reasons that syntax does not convey. Avoid narrating statements.

Review the total source needed to understand a calculation, including helpers
and data structures. A short call hiding machinery is not evidence of clarity.
Review compactness and visible computation independently of correctness and
performance. Correct results alone do not establish style compliance.

Prove zero allocation with instrumentation across all reference phases and
fresh/repeated calls. Verify ordinary and customized imports, caller storage,
foreign layouts, and private-edit rebuild isolation. Preserve scientific
fixtures and historical receipts. Performance and compilation claims require
current, source-bound measurements; source appearance is not a speed result.
