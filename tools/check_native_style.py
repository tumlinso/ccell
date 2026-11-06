#!/usr/bin/env python3
"""Audit shipped ccell and Velocity++ C++ against coding_style.md."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


CPP_SUFFIXES = {
    ".c", ".cc", ".cpp", ".cxx", ".h", ".hh", ".hpp", ".hxx",
    ".cppm", ".ixx", ".mpp", ".cuh", ".inc", ".inl", ".ipp", ".tpp",
}
DEFAULT_SOURCES = [
    "src", "external/velocitypp/src", "examples", "external/velocitypp/examples",
]
FUNDAMENTAL_HEADERS = {"stddef.h", "stdint.h", "float.h"}
INTRINSIC_HEADERS = {
    "arm_neon.h", "avx512bwintrin.h", "avx512dqintrin.h", "avx512fintrin.h",
    "avx512vlintrin.h", "bmiintrin.h", "emmintrin.h", "immintrin.h",
    "pmmintrin.h", "smmintrin.h", "tmmintrin.h", "xmmintrin.h",
}
ANN_DATA_SCALARS = {
    "byte", "size_t", "ptrdiff_t", "max_align_t", "int8_t", "uint8_t",
    "int16_t", "uint16_t", "int32_t", "uint32_t", "int64_t", "uint64_t",
    "intmax_t", "uintmax_t", "uintptr_t", "intptr_t",
}
STD_NAME_RE = re.compile(r"\bstd\s*::\s*([A-Za-z_]\w*)")
STD_NAMESPACE_RE = re.compile(r"\busing\s+namespace\s+std\s*;")
STD_USING_RE = re.compile(r"\busing\s+std\s*::\s*([A-Za-z_]\w*)")
ALLOC_RE = re.compile(
    r"\b(?:new|delete)\b|\b(?:malloc|calloc|realloc|reallocarray|free|"
    r"aligned_alloc|posix_memalign|valloc|strdup|strndup|alloca)\s*\("
)
INCLUDE_RE = re.compile(r"^\s*#\s*include\s*([<\"])([^>\"]+)[>\"]")
PREPROCESSOR_RE = re.compile(r"^\s*#")
STD_NAMESPACE_ALIAS_RE = re.compile(r"\bnamespace\s+[A-Za-z_]\w*\s*=\s*std\s*;")
STD_IMPORT_RE = re.compile(r"\b(?:export\s+)?import\s+(?:std(?:\.[A-Za-z_]\w*)?|<[^>]+>)\s*;")


def mask_non_code(text: str) -> str:
    """Replace comments and literals with spaces while preserving newlines."""
    out = list(text)
    n = len(text)
    i = 0

    def blank(start: int, end: int) -> None:
        for at in range(start, end):
            if out[at] != "\n":
                out[at] = " "

    while i < n:
        if text.startswith("//", i):
            end = text.find("\n", i + 2)
            if end < 0:
                end = n
            blank(i, end)
            i = end
            continue
        if text.startswith("/*", i):
            end = text.find("*/", i + 2)
            end = n if end < 0 else end + 2
            blank(i, end)
            i = end
            continue

        # Raw string literal: R"tag(contents)tag". Prefix forms u8R/uR/UR/LR
        # reach this R through the normal character scan; only the R matters.
        raw = re.match(r'R"([^ ()\\\t\r\n]{0,16})\(', text[i:])
        if raw:
            delimiter = raw.group(1)
            start = i
            content = i + raw.end()
            terminator = ")" + delimiter + '"'
            end = text.find(terminator, content)
            end = n if end < 0 else end + len(terminator)
            blank(start, end)
            i = end
            continue

        # Ordinary and prefixed string/character literals. Prefix letters are
        # left as code identifiers; the quote itself begins the masked region.
        if text[i] in "\"'":
            quote = text[i]
            start = i
            i += 1
            while i < n:
                if text[i] == "\\":
                    i += 2
                    continue
                if text[i] == quote:
                    i += 1
                    break
                i += 1
            blank(start, min(i, n))
            continue
        i += 1
    return "".join(out)


def source_files(root: Path, roots: list[str]) -> list[Path]:
    found: set[Path] = set()
    for rel in roots:
        directory = root / rel
        if directory.is_dir():
            found.update(
                p.resolve() for p in directory.rglob("*")
                if p.is_file() and p.suffix in CPP_SUFFIXES
            )
    return sorted(found)


def under(path: Path, directory: Path) -> bool:
    try:
        path.relative_to(directory)
        return True
    except ValueError:
        return False


def local_include(path: Path, spelling: str, root: Path) -> Path | None:
    for candidate in ((path.parent / spelling).resolve(), (root / spelling).resolve()):
        if under(candidate, root) and candidate.is_file():
            return candidate
    return None


def include_lines(raw: str, masked: str) -> list[tuple[int, str, str]]:
    """Return source include directives as (line, delimiter, path)."""
    raw_lines = raw.splitlines()
    code_lines = masked.splitlines()
    found: list[tuple[int, str, str]] = []
    i = 0
    while i < min(len(raw_lines), len(code_lines)):
        code = code_lines[i]
        raw_line = raw_lines[i]
        if not PREPROCESSOR_RE.match(code):
            i += 1
            continue

        first = i + 1
        logical_raw = raw_line
        logical_code = code
        while logical_raw.rstrip().endswith("\\") and i + 1 < len(raw_lines):
            i += 1
            logical_raw = logical_raw.rstrip()[:-1] + raw_lines[i].lstrip()
            logical_code = logical_code.rstrip()[:-1] + code_lines[i].lstrip()

        if re.match(r"^\s*#\s*include\b", logical_code):
            match = INCLUDE_RE.match(logical_raw)
            if match and not logical_code[match.end():].strip():
                found.append((first, match.group(1), match.group(2).strip()))
            else:
                found.append((first, "?", "<macro or malformed include>"))
        i += 1
    return found


def audit(root: Path, roots: list[str]) -> tuple[list[str], int, int]:
    root = root.resolve()
    authored_root = root / "external/anndatapp"
    reference_roots = (root / "src", root / "external/velocitypp/src")
    errors: list[str] = []
    initial = source_files(root, roots)
    visited: set[Path] = set()

    def scan(path: Path, inherited_allocation_rule: bool) -> None:
        path = path.resolve()
        if under(path, authored_root.resolve()):
            return  # AnnData++ is the explicit imported implementation exception.
        visit_key = (path, inherited_allocation_rule)
        if visit_key in visited or not path.is_file():
            return
        visited.add(visit_key)
        try:
            raw = path.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            return
        code = mask_non_code(raw)
        try:
            rel = path.relative_to(root)
        except ValueError:
            return

        allocation_rule = inherited_allocation_rule or any(
            under(path, directory.resolve()) for directory in reference_roots
        )

        for line, delimiter, include in include_lines(raw, code):
            if delimiter == "<":
                if include not in FUNDAMENTAL_HEADERS and Path(include).name not in INTRINSIC_HEADERS:
                    errors.append(f"{rel}:{line}: system include <{include}> is outside the minimal header set")
                continue
            if delimiter == "?":
                errors.append(f"{rel}:{line}: include path must be explicit so local dependencies can be audited")
                continue
            nested = local_include(path, include, root)
            if nested:
                scan(nested, allocation_rule)
            else:
                errors.append(f"{rel}:{line}: quoted include \"{include}\" is unresolved in authored sources")

        for match in STD_NAME_RE.finditer(code):
            name = match.group(1)
            if name not in ANN_DATA_SCALARS:
                line = code.count("\n", 0, match.start()) + 1
                errors.append(f"{rel}:{line}: standard-library name std::{name} is outside the scalar descriptor exception")
        for match in STD_NAMESPACE_RE.finditer(code):
            line = code.count("\n", 0, match.start()) + 1
            errors.append(f"{rel}:{line}: importing the std namespace is prohibited")
        for match in STD_NAMESPACE_ALIAS_RE.finditer(code):
            line = code.count("\n", 0, match.start()) + 1
            errors.append(f"{rel}:{line}: aliasing the std namespace is prohibited")
        for match in STD_IMPORT_RE.finditer(code):
            line = code.count("\n", 0, match.start()) + 1
            errors.append(f"{rel}:{line}: importing the standard library module/header is prohibited")
        for match in STD_USING_RE.finditer(code):
            name = match.group(1)
            if name in ANN_DATA_SCALARS:
                continue
            line = code.count("\n", 0, match.start()) + 1
            errors.append(f"{rel}:{line}: importing std::{name} is outside the scalar descriptor exception")

        if allocation_rule:
            for match in ALLOC_RE.finditer(code):
                token = match.group(0)
                if token == "delete":
                    # `= delete` is a special-member declaration, not an action.
                    line_start = code.rfind("\n", 0, match.start()) + 1
                    line_end = code.find("\n", match.start())
                    line_end = len(code) if line_end < 0 else line_end
                    if re.search(r"=\s*$", code[line_start:match.start()]):
                        continue
                line = code.count("\n", 0, match.start()) + 1
                errors.append(f"{rel}:{line}: dynamic allocation/deallocation in the standard reference path")

    for path in initial:
        scan(path, False)
    return sorted(set(errors)), len(initial), len({path for path, _ in visited})


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument(
        "--source", action="append", dest="sources",
        help="authored shipped C++ root to audit (repeatable; defaults to production roots)",
    )
    args = parser.parse_args()
    root = args.root.resolve()
    roots = args.sources or DEFAULT_SOURCES
    errors, source_count, scanned_count = audit(root, roots)
    if errors:
        print("\n".join(errors), file=sys.stderr)
        print(f"ccell native style failed: {len(errors)} finding(s), {scanned_count} files scanned", file=sys.stderr)
        return 1
    print(
        f"ccell native style passed: {source_count} shipped C++ files; "
        f"{scanned_count} files including authored local headers"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
