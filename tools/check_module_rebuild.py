#!/usr/bin/env python3
"""Check ccell module rebuild boundaries in a disposable full-source copy.

The active checkout is read-only to this tool. All copied sources, build trees,
logs, and the JSON report are written below the checkout's build directory.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class Stamp:
    sha256: str
    mtime_ns: int
    size: int


class CheckFailure(RuntimeError):
    pass


IGNORED_SOURCE_NAMES = {
    ".git", ".ccell-local", ".codex", ".agents", "build", "CMakeFiles",
    "__pycache__", ".pytest_cache", ".mypy_cache", ".ruff_cache", ".cache",
    ".venv", "venv", "virtualenv", ".virtualenv", ".tox", ".nox",
    "cmake-build-debug", "cmake-build-release", "Testing", "dist", "htmlcov",
}


def ignored_source_name(name: str) -> bool:
    return (
        name in IGNORED_SOURCE_NAMES
        or name.startswith("cmake-build-")
        or name.endswith(".egg-info")
    )


def within(path: Path, parent: Path) -> bool:
    try:
        path.resolve().relative_to(parent.resolve())
        return True
    except ValueError:
        return False


def parse_args() -> argparse.Namespace:
    source = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=source,
                        help="current checkout to copy (defaults to this script's repository)")
    parser.add_argument("--artifacts-dir", type=Path,
                        help="artifact/work directory under SOURCE/build (default: timestamped directory)")
    parser.add_argument("--output", type=Path,
                        help="JSON report path under SOURCE/build (default: ARTIFACTS_DIR/report.json)")
    parser.add_argument("--compiler", default="clang++-18")
    parser.add_argument("--scanner", default="clang-scan-deps-18")
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--timeout-seconds", type=float, default=1800,
                        help="per-subprocess timeout; must exceed 30 seconds (default: 1800)")
    args = parser.parse_args()
    args.source_dir = args.source_dir.resolve()
    build_root = (args.source_dir / "build").resolve()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
    args.artifacts_dir = (args.artifacts_dir or build_root / "artifacts" / "module-rebuild" / stamp).resolve()
    args.output = (args.output or args.artifacts_dir / "report.json").resolve()
    if args.jobs < 1 or args.jobs > 2:
        parser.error("--jobs must be 1 or 2; this check caps build concurrency at two")
    if not math.isfinite(args.timeout_seconds) or args.timeout_seconds <= 30:
        parser.error("--timeout-seconds must be greater than 30")
    if not within(args.artifacts_dir, build_root) or not within(args.output, build_root):
        parser.error("--artifacts-dir and --output must be under SOURCE/build")
    if args.artifacts_dir == build_root or args.output == build_root:
        parser.error("--artifacts-dir and --output must be below SOURCE/build")
    if args.artifacts_dir.exists():
        parser.error(f"refusing to reuse existing artifacts directory: {args.artifacts_dir}")
    if args.output.exists():
        parser.error(f"refusing to overwrite existing report: {args.output}")
    if not args.source_dir.is_dir():
        parser.error(f"source directory does not exist: {args.source_dir}")
    return args


def copy_source(source: Path, destination: Path) -> None:
    """Copy all current source, including untracked edits, without Git/build output."""
    if destination.exists():
        raise CheckFailure(f"refusing to reuse existing disposable source tree: {destination}")
    links: list[Path] = []
    for directory, child_dirs, files in os.walk(source, followlinks=False):
        base = Path(directory)
        child_dirs[:] = [name for name in child_dirs if not ignored_source_name(name)]
        links.extend(base / name for name in child_dirs if (base / name).is_symlink())
        links.extend(
            base / name for name in files
            if not ignored_source_name(name) and (base / name).is_symlink()
        )
        child_dirs[:] = [name for name in child_dirs if not (base / name).is_symlink()]
    if links:
        details = ", ".join(str(p.relative_to(source)) for p in links[:8])
        raise CheckFailure(
            "source copy contains symlinks; resolve or copy their targets explicitly before running: "
            + details
        )

    def ignore(_directory: str, names: list[str]) -> set[str]:
        return {name for name in names if ignored_source_name(name)}

    shutil.copytree(source, destination, ignore=ignore, symlinks=False)


def write_source_manifest(source: Path, output: Path) -> dict[str, Any]:
    entries = []
    total_bytes = 0
    for directory, child_dirs, files in os.walk(source, followlinks=False):
        base = Path(directory)
        child_dirs[:] = [name for name in child_dirs if not ignored_source_name(name)]
        for name in files:
            if ignored_source_name(name):
                continue
            path = base / name
            if not path.is_file():
                continue
            size = path.stat().st_size
            total_bytes += size
            entries.append({
                "path": path.relative_to(source).as_posix(),
                "size": size,
                "sha256": sha256(path),
            })
    entries.sort(key=lambda item: item["path"])
    encoded = json.dumps(entries, sort_keys=True, separators=(",", ":")).encode()
    output.write_bytes(encoded + b"\n")
    return {
        "path": str(output), "sha256": hashlib.sha256(encoded).hexdigest(),
        "file_count": len(entries), "total_bytes": total_bytes,
    }


def run(command: list[str], cwd: Path, log: Path, timeout_seconds: float,
        env: dict[str, str] | None = None,
        expected: set[int] | frozenset[int] = frozenset({0})) -> tuple[int, float, str]:
    start = time.monotonic()
    process = subprocess.Popen(
        command, cwd=cwd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        env=env, start_new_session=(os.name == "posix"),
    )
    try:
        stdout, stderr = process.communicate(timeout=timeout_seconds)
    except subprocess.TimeoutExpired:
        if os.name == "posix":
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            try:
                stdout, stderr = process.communicate(timeout=5)
            except subprocess.TimeoutExpired:
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                stdout, stderr = process.communicate()
        else:
            process.kill()
            stdout, stderr = process.communicate()
        elapsed = time.monotonic() - start
        output = (
            f"$ {' '.join(command)}\n\n[timeout after {timeout_seconds:g} seconds]\n"
            f"[stdout]\n{stdout or ''}\n[stderr]\n{stderr or ''}"
        )
        log.parent.mkdir(parents=True, exist_ok=True)
        log.write_text(output, encoding="utf-8")
        raise CheckFailure(
            f"command exceeded {timeout_seconds:g} seconds: {' '.join(command)}\nlog: {log}"
        )
    elapsed = time.monotonic() - start
    output = f"$ {' '.join(command)}\n\n[stdout]\n{stdout}\n[stderr]\n{stderr}"
    log.parent.mkdir(parents=True, exist_ok=True)
    log.write_text(output, encoding="utf-8")
    if process.returncode not in expected:
        raise CheckFailure(
            f"command failed with exit {process.returncode}: {' '.join(command)}\n"
            f"log: {log}\n{stdout[-5000:]}\n{stderr[-5000:]}"
        )
    return process.returncode, elapsed, stdout + stderr


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def snapshots(build: Path, suffixes: set[str]) -> dict[str, Stamp]:
    result: dict[str, Stamp] = {}
    if not build.exists():
        return result
    for path in build.rglob("*"):
        if not path.is_file() or path.suffix not in suffixes:
            continue
        try:
            rel = path.relative_to(build).as_posix()
            stat = path.stat()
            result[rel] = Stamp(sha256(path), stat.st_mtime_ns, stat.st_size)
        except OSError:
            continue
    return result


def changed(before: dict[str, Stamp], after: dict[str, Stamp]) -> list[str]:
    keys = before.keys() | after.keys()
    return sorted(key for key in keys if before.get(key) != after.get(key))


def stamp_manifest(items: dict[str, Stamp]) -> dict[str, dict[str, Any]]:
    return {path: asdict(items[path]) for path in sorted(items)}


def manifest_hash(items: dict[str, Stamp]) -> str:
    encoded = json.dumps(stamp_manifest(items), sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def consumer_object(path: str) -> bool:
    if not path.endswith(".o"):
        return False
    target = re.search(r"(?:^|/)CMakeFiles/([^/]+)\.dir/", path)
    if not target:
        return False
    name = target.group(1)
    if "installed/consumer/CMakeFiles/consumer_" in path:
        return True
    module_providers = {"anndatapp_core", "anndatapp", "ccell_preprocess", "velocitypp", "ccell"}
    return name not in module_providers


def module_interface_object(path: str) -> bool:
    return path.endswith(".cppm.o") or path.endswith(".ixx.o")


def public_bmi(path: str) -> bool:
    return path.endswith(".pcm") or path.endswith(".ifc")


def build_config(source: Path, build: Path, config: str, compiler: str,
                 scanner: str, jobs: int, timeout_seconds: float,
                 logs: Path) -> dict[str, Any]:
    compiler_path = shutil.which(compiler)
    scanner_path = shutil.which(scanner)
    if not compiler_path:
        raise CheckFailure(f"compiler not found: {compiler}")
    if not scanner_path:
        raise CheckFailure(f"matching module dependency scanner not found: {scanner}")
    build.mkdir(parents=True, exist_ok=True)
    configure = [
        "cmake", "-S", str(source), "-B", str(build), "-G", "Ninja",
        f"-DCMAKE_CXX_COMPILER={compiler_path}",
        f"-DCMAKE_CXX_COMPILER_CLANG_SCAN_DEPS={scanner_path}",
        f"-DCMAKE_BUILD_TYPE={config}",
        "-DANNDATAPP_ENABLE_HDF5=OFF", "-DANNDATAPP_ENABLE_ZARR=OFF",
        "-DCCELL_BUILD_TESTS=ON",
    ]
    _, configure_s, _ = run(configure, source, logs / f"{config}-configure.log", timeout_seconds)
    build_cmd = ["cmake", "--build", str(build), "--parallel", str(jobs)]
    _, build_s, _ = run(build_cmd, source, logs / f"{config}-build.log", timeout_seconds)
    _, test_s, _ = run(
        ["ctest", "--test-dir", str(build), "--output-on-failure"],
        source, logs / f"{config}-ctest.log", timeout_seconds,
    )
    return {
        "build_dir": str(build), "configure_seconds": configure_s,
        "build_seconds": build_s, "ctest_seconds": test_s,
        "tests": "passed",
    }


def append_probe_comment(path: Path, label: str) -> tuple[str, str]:
    before_hash = sha256(path)
    with path.open("a", encoding="utf-8") as stream:
        stream.write(f"\n// ccell module rebuild probe: {label}\n")
    os.utime(path, None)
    return before_hash, sha256(path)


def private_edit(build: Path, source_file: Path, label: str, logs: Path,
                 jobs: int, timeout_seconds: float) -> dict[str, Any]:
    before_all = snapshots(build, {".pcm", ".ifc", ".o"})
    before_bmi = {p: s for p, s in before_all.items() if public_bmi(p)}
    before_consumer = {p: s for p, s in before_all.items() if consumer_object(p)}
    before_interface_objects = {p: s for p, s in before_all.items() if module_interface_object(p)}
    source_hash_before, source_hash_after = append_probe_comment(source_file, label)
    _, elapsed, output = run(
        ["cmake", "--build", str(build), "--parallel", str(jobs)],
        build.parent, logs, timeout_seconds,
    )
    after_all = snapshots(build, {".pcm", ".ifc", ".o"})
    after_bmi = {p: s for p, s in after_all.items() if public_bmi(p)}
    after_consumer = {p: s for p, s in after_all.items() if consumer_object(p)}
    after_interface_objects = {p: s for p, s in after_all.items() if module_interface_object(p)}
    bmi_changed = changed(before_bmi, after_bmi)
    consumer_changed = changed(before_consumer, after_consumer)
    interface_changed = changed(before_interface_objects, after_interface_objects)
    object_changed = [p for p in changed(before_all, after_all) if p.endswith(".o")]
    if bmi_changed:
        raise CheckFailure(f"private edit {label} changed public BMI artifacts: {bmi_changed}")
    if consumer_changed:
        raise CheckFailure(f"private edit {label} rebuilt consumer objects: {consumer_changed}")
    if interface_changed:
        raise CheckFailure(f"private edit {label} rebuilt module interface objects: {interface_changed}")
    if not object_changed:
        raise CheckFailure(f"private edit {label} did not rebuild any implementation object")
    expected_object = source_file.name + ".o"
    if not any(Path(path).name == expected_object for path in object_changed):
        raise CheckFailure(
            f"private edit {label} did not rebuild the object for {source_file.name}: {object_changed}"
        )
    return {
        "source": str(source_file), "seconds": elapsed,
        "source_sha256_before": source_hash_before,
        "source_sha256_after_probe_edit": source_hash_after,
        "changed_objects": object_changed, "bmi_unchanged": True,
        "consumer_objects_unchanged": True, "module_interface_objects_unchanged": True,
        "public_bmis_before": stamp_manifest(before_bmi),
        "public_bmis_after": stamp_manifest(after_bmi),
        "consumer_objects_before": stamp_manifest(before_consumer),
        "consumer_objects_after": stamp_manifest(after_consumer),
        "module_interface_objects_before": stamp_manifest(before_interface_objects),
        "module_interface_objects_after": stamp_manifest(after_interface_objects),
        "build_output_tail": output[-1800:],
    }


def consumer_edit(build: Path, source_file: Path, logs: Path, jobs: int,
                  timeout_seconds: float) -> dict[str, Any]:
    before_all = snapshots(build, {".pcm", ".ifc", ".o"})
    before_bmi = {p: s for p, s in before_all.items() if public_bmi(p)}
    before_objects = {p: s for p, s in before_all.items() if p.endswith(".o")}
    before_interface_objects = {p: s for p, s in before_all.items() if module_interface_object(p)}
    source_hash_before, source_hash_after = append_probe_comment(source_file, "consumer source edit")
    _, build_seconds, output = run(
        ["cmake", "--build", str(build), "--parallel", str(jobs)],
        build.parent, logs / "consumer-source-build.log", timeout_seconds,
    )
    _, ctest_seconds, ctest_output = run(
        ["ctest", "--test-dir", str(build), "-R", "ccell-installed-consumer", "--output-on-failure"],
        build.parent, logs / "consumer-source-installed-ctest.log", timeout_seconds,
    )
    after_all = snapshots(build, {".pcm", ".ifc", ".o"})
    after_bmi = {p: s for p, s in after_all.items() if public_bmi(p)}
    after_objects = {p: s for p, s in after_all.items() if p.endswith(".o")}
    after_interface_objects = {p: s for p, s in after_all.items() if module_interface_object(p)}
    bmi_changed = changed(before_bmi, after_bmi)
    objects_changed = changed(before_objects, after_objects)
    interface_changed = changed(before_interface_objects, after_interface_objects)
    unexpected = [p for p in objects_changed if not consumer_object(p)]
    if bmi_changed:
        raise CheckFailure(f"editing a consumer source changed module BMIs: {bmi_changed}")
    if interface_changed:
        raise CheckFailure(f"editing a consumer source rebuilt module interface objects: {interface_changed}")
    if unexpected:
        raise CheckFailure(f"editing a consumer source rebuilt non-consumer objects: {unexpected}")
    if not objects_changed:
        raise CheckFailure("consumer source edit did not rebuild a consumer object")
    return {
        "source": str(source_file), "changed_consumer_objects": objects_changed,
        "source_sha256_before": source_hash_before,
        "source_sha256_after_probe_edit": source_hash_after,
        "module_bmis_unchanged": True, "other_objects_unchanged": True,
        "module_interface_objects_unchanged": True,
        "public_bmis_before": stamp_manifest(before_bmi),
        "public_bmis_after": stamp_manifest(after_bmi),
        "all_object_files_before": stamp_manifest(before_objects),
        "all_object_files_after": stamp_manifest(after_objects),
        "module_interface_objects_before": stamp_manifest(before_interface_objects),
        "module_interface_objects_after": stamp_manifest(after_interface_objects),
        "build_seconds": build_seconds, "installed_consumer_seconds": ctest_seconds,
        "build_output_tail": output[-1800:], "ctest_output_tail": ctest_output[-1800:],
    }


def specialization_probe(artifacts: Path, compiler: str, scanner: str,
                         jobs: int, timeout_seconds: float,
                         logs: Path) -> dict[str, Any]:
    source = artifacts / "specialization-probe-src"
    build = artifacts / "specialization-probe-build"
    source.mkdir(parents=True)
    (source / "CMakeLists.txt").write_text(
        "cmake_minimum_required(VERSION 3.28)\n"
        "project(ccell_rebuild_specialization LANGUAGES CXX)\n"
        "set(CMAKE_CXX_EXTENSIONS OFF)\n"
        "add_library(rebuild_probe)\n"
        "target_compile_features(rebuild_probe PUBLIC cxx_std_20)\n"
        "target_sources(rebuild_probe PUBLIC FILE_SET CXX_MODULES FILES probe.cppm)\n"
        "add_executable(specialization_consumer main.cpp)\n"
        "target_link_libraries(specialization_consumer PRIVATE rebuild_probe)\n",
        encoding="utf-8",
    )
    module = source / "probe.cppm"
    module.write_text(
        "export module ccell_rebuild_specialization;\n"
        "export template<int N> constexpr int value() { return N; }\n"
        "export template<> constexpr int value<7>() { return 7; }\n",
        encoding="utf-8",
    )
    (source / "main.cpp").write_text(
        "import ccell_rebuild_specialization;\n"
        "int main() { return value<7>(); }\n",
        encoding="utf-8",
    )
    compiler_path = shutil.which(compiler)
    scanner_path = shutil.which(scanner)
    if not compiler_path or not scanner_path:
        raise CheckFailure("compiler or module dependency scanner disappeared before specialization probe")
    _, configure_seconds, _ = run(
        ["cmake", "-S", str(source), "-B", str(build), "-G", "Ninja",
         f"-DCMAKE_CXX_COMPILER={compiler_path}",
         f"-DCMAKE_CXX_COMPILER_CLANG_SCAN_DEPS={scanner_path}",
         "-DCMAKE_BUILD_TYPE=Debug"],
        source, logs / "specialization-configure.log", timeout_seconds,
    )
    command = ["cmake", "--build", str(build), "--parallel", str(jobs)]
    _, first_build_seconds, _ = run(
        command, source, logs / "specialization-initial-build.log", timeout_seconds,
    )
    exe = build / "specialization_consumer"
    initial_result, _, _ = run(
        [str(exe)], build, logs / "specialization-initial-run.log", timeout_seconds,
        expected={7},
    )
    before = snapshots(build, {".pcm", ".ifc", ".o"})
    original = module.read_text(encoding="utf-8")
    updated = original.replace("return 7;", "return 8;")
    if updated == original:
        raise CheckFailure("specialization probe source did not contain its expected edit point")
    module.write_text(updated, encoding="utf-8")
    os.utime(module, None)
    _, changed_build_seconds, output = run(
        command, source, logs / "specialization-changed-build.log", timeout_seconds,
    )
    after = snapshots(build, {".pcm", ".ifc", ".o"})
    changed_bmis = [p for p in changed(before, after) if public_bmi(p)]
    changed_consumers = [p for p in changed(before, after) if p.endswith(".o") and "main.cpp" in p]
    if not changed_bmis:
        raise CheckFailure("specialization edit did not rebuild its module BMI")
    if not changed_consumers:
        raise CheckFailure("specialization edit did not rebuild the importing consumer object")
    final_result, _, _ = run(
        [str(exe)], build, logs / "specialization-changed-run.log", timeout_seconds,
        expected={8},
    )
    return {
        "configure_seconds": configure_seconds,
        "initial_build_seconds": first_build_seconds,
        "changed_build_seconds": changed_build_seconds,
        "initial_result": initial_result,
        "changed_result": final_result,
        "changed_module_bmis": changed_bmis,
        "changed_consumer_objects": changed_consumers,
        "build_output_tail": output[-1800:],
    }


def choose_source(paths: list[Path], description: str) -> Path:
    if not paths:
        raise CheckFailure(f"no private {description} implementation .cpp file found in this source layout")
    return sorted(paths)[0]


def main() -> int:
    args = parse_args()
    args.artifacts_dir.mkdir(parents=True, exist_ok=False)
    source_copy = args.artifacts_dir / "source"
    report: dict[str, Any] = {
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "source_checkout": str(args.source_dir),
        "source_copy": str(source_copy),
        "compiler": args.compiler,
        "scanner": args.scanner,
        "jobs": args.jobs,
        "timeout_seconds": args.timeout_seconds,
        "artifacts_dir": str(args.artifacts_dir),
        "checks": {},
    }
    try:
        print(f"Copying current sources to {source_copy}", flush=True)
        copy_source(args.source_dir, source_copy)
        report["source_snapshot"] = write_source_manifest(
            source_copy, args.artifacts_dir / "source-snapshot-initial.json",
        )
        builds = args.artifacts_dir / "builds"
        logs = args.artifacts_dir / "logs"
        report["configurations"] = {}
        for config in ("Debug", "Release"):
            build = builds / config.lower()
            print(f"Configuring, building, and testing {config}", flush=True)
            report["configurations"][config.lower()] = build_config(
                source_copy, build, config, args.compiler, args.scanner, args.jobs,
                args.timeout_seconds, logs,
            )

        debug = builds / "debug"
        print("Checking no-op Debug rebuild", flush=True)
        initial = snapshots(debug, {".pcm", ".ifc", ".o"})
        _, noop_seconds, noop_output = run(
            ["cmake", "--build", str(debug), "--parallel", str(args.jobs)],
            source_copy, logs / "debug-noop-build.log", args.timeout_seconds,
        )
        after_noop = snapshots(debug, {".pcm", ".ifc", ".o"})
        noop_changes = changed(initial, after_noop)
        if noop_changes:
            raise CheckFailure(f"no-op Debug build modified compile artifacts: {noop_changes}")
        if "no work to do" not in noop_output.lower():
            raise CheckFailure("Ninja did not report a no-op build; inspect debug-noop-build.log")
        report["checks"]["debug_noop_build"] = {
            "seconds": noop_seconds, "compile_artifacts_unchanged": True,
            "compile_artifact_count": len(initial),
            "compile_artifact_manifest_sha256": manifest_hash(initial),
            "output": noop_output[-800:],
        }

        umbrella = source_copy / "tests/umbrella.cpp"
        if not umbrella.is_file():
            raise CheckFailure("tests/umbrella.cpp is required for the consumer-only edit probe")
        print("Checking consumer-only source rebuild", flush=True)
        report["checks"]["consumer_source_edit"] = consumer_edit(
            debug, umbrella, logs, args.jobs, args.timeout_seconds,
        )

        preprocess_sources = [
            p for p in (source_copy / "src/preprocess").rglob("*.cpp")
            if p.is_file() and "test" not in p.parts and "tests" not in p.parts
        ]
        velocity_sources = [
            p for p in (source_copy / "external/velocitypp/src").rglob("*.cpp")
            if p.is_file() and "test" not in p.parts and "tests" not in p.parts
        ]
        preprocess = choose_source(preprocess_sources, "preprocess")
        velocity = choose_source(velocity_sources, "Velocity++")
        print(f"Checking private preprocess edit: {preprocess.relative_to(source_copy)}", flush=True)
        report["checks"]["private_preprocess_edit"] = private_edit(
            debug, preprocess, "preprocess", logs / "private-preprocess-build.log",
            args.jobs, args.timeout_seconds,
        )
        print(f"Checking private Velocity++ edit: {velocity.relative_to(source_copy)}", flush=True)
        report["checks"]["private_velocity_edit"] = private_edit(
            debug, velocity, "velocity", logs / "private-velocity-build.log",
            args.jobs, args.timeout_seconds,
        )
        print("Checking compile-time specialization rebuild", flush=True)
        report["checks"]["specialization_probe"] = specialization_probe(
            args.artifacts_dir, args.compiler, args.scanner, args.jobs,
            args.timeout_seconds, logs,
        )
        report["source_snapshot_after_probe_edits"] = write_source_manifest(
            source_copy, args.artifacts_dir / "source-snapshot-after-edits.json",
        )
        report["status"] = "passed"
    except Exception as exc:
        report["status"] = "failed"
        report["error"] = f"{type(exc).__name__}: {exc}"
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(report["error"], file=sys.stderr)
        print(f"partial report: {args.output}", file=sys.stderr)
        return 1

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"module rebuild contract passed; report: {args.output}")
    print(f"disposable sources and build artifacts: {args.artifacts_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
