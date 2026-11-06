# Working on ccell

ccell owns a C++20 reference biology package. AnnData++ and Velocity++ are
standalone dependencies included as Git submodules; preprocessing is local
ccell source under src/preprocess. Dependencies live in external/. Cellerator is a downstream consumer and the home of FP16 and
specialized device implementations. Do not add a Cellerator dependency to
this reference package.

Read [coding_style.md](coding_style.md) before numerical design, implementation,
or review. It is binding authority: do not weaken, replace, move, or remove it,
or its protection here, without explicit user authorization. The standard
reference allocates nothing; callers supply all outputs, preparation, and scratch.
Consumer replacements may allocate explicitly. Shipped source and examples use
minimal scalar headers and compiler builtins, without STL machinery.

Public imports are `ccell`, `ccell.anndata`, `ccell.preprocess`, and
`ccell.velocity`. Re-export the existing component modules and preserve
their type identity. AnnData++ owns arrays, views, layouts, and lifetimes;
keep numerical operations native and memory/allocation behavior explicit.
Keep implementation-only changes out of public module interfaces.

Use the debug preset with Clang 18, matching clang-scan-deps, CMake >=3.28,
and Ninja: configure, build, then run CTest. It checks composition through
the umbrella and individual module imports and through an installed package.
The same root build owns preprocessing's numerical tests and packaging; do
not recreate standalone preprocess build or package metadata. Its license
and third-party notices live at the repository root. Velocity++ numerical
changes also require its standalone checks. Historical qualification records
are not current-source evidence.

Preserve concurrent and unknown work, staged changes, submodule histories,
scientific fixtures, and historical receipts. Do not reset or clean work,
delete remote repositories, or rewrite workflow authority. Use configured
Codex subagents for useful independent work, with explicit file ownership;
keep intent, integration, and acceptance with the root.
