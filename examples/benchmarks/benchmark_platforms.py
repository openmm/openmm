#!/usr/bin/env python3
# Copyright (c) 2026 Chun-Chi Hung.
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

"""Record benchmark.py results for selected independently built platform variants.

This driver uses only the standard library; workers require OpenMM Python bindings
built from this checkout. No packages are installed by the driver. Run --dry-run
to inspect the plan without building, importing OpenMM, or accessing a GPU.
"""

import argparse
import ctypes
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import runpy
import shutil
import signal
import statistics
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
BENCHMARK = ROOT / "examples/benchmarks/benchmark.py"
TESTS = ("gbsa", "rf", "pme", "apoa1rf", "apoa1pme", "apoa1ljpme")
VARIANTS = ("opencl", "metal_off", "metal_on")


def digest(path):
    """Hash inputs and actual binaries, not merely CMake option labels."""
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            value.update(block)
    return value.hexdigest()


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+"\n")


def cache_values(build):
    result = {}
    for line in (build / "CMakeCache.txt").read_text().splitlines():
        if not line or line.startswith(("#", "//")) or "=" not in line:
            continue
        name, value = line.split("=", 1)
        result[name.split(":", 1)[0]] = value
    return result


def configurations(args):
    """Reject mislabeled comparison endpoints before launching any benchmark."""
    configs = []
    for label, platform, directory in (
            ("opencl", "OpenCL", args.opencl_build),
            ("metal_off", "Metal", args.metal_off_build),
            ("metal_on", "Metal", args.metal_on_build)):
        if label not in args.variants:
            continue
        build = Path(directory).expanduser().resolve()
        cache = cache_values(build)
        if cache.get("CMAKE_BUILD_TYPE") != "Release":
            raise ValueError(f"{label}: expected a Release build: {build}")
        if Path(cache.get("CMAKE_HOME_DIRECTORY", "")).resolve() != ROOT:
            raise ValueError(f"{label}: build does not belong to this checkout")
        if cache.get(f"OPENMM_BUILD_{platform.upper()}_LIB") != "ON":
            raise ValueError(f"{label}: platform library is not enabled")
        flags = {key: value for key, value in cache.items()
                 if key.startswith("OPENMM_METAL_") and value in ("ON", "OFF")}
        if platform == "Metal":
            fast = {key: value for key, value in flags.items() if key.startswith("OPENMM_METAL_FAST_")}
            if not fast:
                raise ValueError(f"{label}: no Metal optimization switches found")
            expected = "ON" if label == "metal_on" else "OFF"
            bad = [key for key, value in fast.items() if value != expected]
            if flags.get("OPENMM_METAL_NATIVE_FLOAT_ATOMICS") != expected:
                bad.append("OPENMM_METAL_NATIVE_FLOAT_ATOMICS")
            # Batching is enabled by turning per-operation commit OFF.
            commit = "OFF" if label == "metal_on" else "ON"
            if flags.get("OPENMM_METAL_RECORD_AND_COMMIT") != commit:
                bad.append("OPENMM_METAL_RECORD_AND_COMMIT")
            for experimental in ("OPENMM_METAL_EXPERIMENTAL_MATRIX_SCREEN", "OPENMM_METAL_EXPERIMENTAL_NONBONDED_HYBRID"):
                if flags.get(experimental, "OFF") != "OFF":
                    bad.append(experimental)
            # Launch-geometry and language-target experiments need their own labels, even when all
            # ordinary fast paths match the standard ON/OFF endpoint.
            bad.extend(key for key, value in flags.items()
                       if (key.startswith("OPENMM_METAL_TUNE_FORCE_") or
                           key == "OPENMM_METAL_TUNE_LANGUAGE_VERSION") and value != "OFF")
            if bad:
                raise ValueError(f"{label}: wrong settings: {', '.join(bad)}")
        configs.append(dict(label=label, platform=platform, build=str(build),
            library=str(build / f"libOpenMM{platform}.dylib"), core=str(build / "libOpenMM.dylib"),
            flags=flags, cmake=cache))
    selected = {config["label"]: config for config in configs}
    if "metal_off" in selected and "metal_on" in selected:
        if {key for key in selected["metal_off"]["flags"] if key.startswith("OPENMM_METAL_FAST_")} != {
                key for key in selected["metal_on"]["flags"] if key.startswith("OPENMM_METAL_FAST_")}:
            raise ValueError("Metal builds expose different fast-path switches; reconfigure both")
    return configs


def run_logged(command, log, environment=None):
    """Stream output and preserve it; interruption stops this exact child group."""
    start = time.monotonic()
    with Path(log).open("w") as output:
        child = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, env=environment, start_new_session=True)
        try:
            for line in child.stdout:
                output.write(line)
                output.flush()
                print(line, end="", flush=True)
            code = child.wait()
        except BaseException:
            if child.poll() is None:
                os.killpg(child.pid, signal.SIGINT)
                try:
                    child.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    os.killpg(child.pid, signal.SIGKILL)
                    child.wait()
            raise
    return code, time.monotonic()-start


def loaded_core():
    """Resolve the actual macOS core library, including pip/conda loader paths."""
    dyld = ctypes.CDLL(None)
    dyld._dyld_image_count.restype = ctypes.c_uint32
    dyld._dyld_get_image_name.argtypes = [ctypes.c_uint32]
    dyld._dyld_get_image_name.restype = ctypes.c_char_p
    images = {Path(os.fsdecode(dyld._dyld_get_image_name(i))).resolve()
              for i in range(dyld._dyld_image_count())}
    cores = [path for path in images if path.name == "libOpenMM.dylib"]
    if len(cores) != 1:
        raise RuntimeError(f"Expected one loaded libOpenMM.dylib, found {cores}")
    return cores[0]


def worker(job_path):
    """Load one platform, fix initialization seeds, then execute the original script."""
    job = json.loads(Path(job_path).read_text())
    try:
        import openmm as mm
    except ImportError as error:
        raise RuntimeError("Install/build the OpenMM Python bindings from this checkout, then "
                           "select that environment with --python. This driver installs nothing.") from error
    package_revision = getattr(mm.version, "git_revision", None)
    if package_revision != job["revision"]:
        raise RuntimeError(f"Python bindings revision {package_revision!r} does not match "
                           f"checkout {job['revision']}; use bindings built from this checkout")
    core = loaded_core()
    if core != Path(job["config"]["core"]).resolve():
        raise RuntimeError(f"Wrong core library loaded: {core}; expected {job['config']['core']}")
    names = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]
    if names != ["Reference"]:
        raise RuntimeError(f"Unexpected auto-loaded platforms: {names}; use an empty plugin directory")
    mm.Platform.loadPluginLibrary(job["config"]["library"])
    if loaded_core() != core:
        raise RuntimeError("Loading the platform changed the core library")
    for key in ("core", "library"):
        if digest(job["config"][key]) != job["config"]["binary_sha256"][key]:
            raise RuntimeError(f"{key} binary changed after the run was configured; stop concurrent builds")
    platform = mm.Platform.getPlatformByName(job["config"]["platform"])
    if "UseCpuPme" in platform.getPropertyNames():
        platform.setPropertyDefaultValue("UseCpuPme", "false")
    # These wrappers change only seed selection, never force or timing code.
    original_init = mm.LangevinMiddleIntegrator.__init__

    def seeded_init(self, *args, **kwargs):
        original_init(self, *args, **kwargs)
        self.setRandomNumberSeed(job["seed"])

    mm.LangevinMiddleIntegrator.__init__ = seeded_init
    original_velocities = mm.Context.setVelocitiesToTemperature

    def seeded_velocities(self, temperature, randomSeed=0):
        return original_velocities(self, temperature, randomSeed or job["seed"])

    mm.Context.setVelocitiesToTemperature = seeded_velocities
    metadata = dict(python=sys.executable, openmm_package=mm.__file__,
        python_extension=mm._openmm.__file__, bindings_revision=package_revision,
        core=str(core), core_sha256=digest(core), plugin=job["config"]["library"],
        plugin_sha256=digest(job["config"]["library"]),
        openmm_version=mm.Platform.getOpenMMVersion(), seed=job["seed"])
    save(Path(job_path).with_name("runtime.json"), metadata)
    if job["test"] is None:
        print(json.dumps(metadata, indent=2))
        return
    os.chdir(BENCHMARK.parent)
    sys.argv = [str(BENCHMARK), "--platform", job["config"]["platform"],
        "--test", job["test"], "--precision", "single", "--ensemble", "NVT",
        "--bond-constraints", "hbonds", "--pme-cutoff", "0.9", "--device", job["device"],
        "--seconds", str(job["seconds"]), "--style", "simple", "--verbose",
        "--outfile", str(Path(job_path).with_name("result.json"))]
    runpy.run_path(str(BENCHMARK), run_name="__main__")


def validate_result(path, job):
    """benchmark.py can swallow an exception and exit zero without a result."""
    rows = json.loads(Path(path).read_text()).get("benchmarks", [])
    if len(rows) != 1:
        raise ValueError(f"Expected exactly one benchmark result, got {len(rows)}")
    row = rows[0]
    for key, expected in (("test", job["test"]), ("platform", job["config"]["platform"]),
                          ("precision", "single"), ("ensemble", "NVT"), ("constraints", "HBonds")):
        if row.get(key) != expected:
            raise ValueError(f"Wrong {key}: {row.get(key)!r}, expected {expected!r}")
    for key in ("steps", "elapsed_time", "ns_per_day"):
        if type(row.get(key)) not in (int, float) or not math.isfinite(row[key]) or row[key] <= 0:
            raise ValueError(f"Invalid {key}: {row.get(key)!r}")
    if type(row["steps"]) is not int:
        raise ValueError("Step count must be an integer")
    cutoff = 2.0 if job["test"] == "gbsa" else 1.0 if job["test"] in ("rf", "apoa1rf") else 0.9
    if row.get("timestep_in_fs") != 4 or row.get("hydrogen_mass") != "1.5" or row.get("cutoff") != cutoff:
        raise ValueError("Unexpected timestep, hydrogen mass, or cutoff")
    expected_speed = row["timestep_in_fs"]*0.0864*row["steps"]/row["elapsed_time"]
    if not math.isclose(row["ns_per_day"], expected_speed, rel_tol=1e-9):
        raise ValueError("ns/day is inconsistent with timestep, elapsed time, and step count")
    props = row["platform_properties"]
    if props.get("Precision") != "single" or props.get("DeviceIndex") != job["device"]:
        raise ValueError(f"Wrong device/precision properties: {props}")
    if str(props.get("UseCpuPme", "false")).lower() != "false":
        raise ValueError("CPU PME was enabled")
    if "Apple" not in props.get("DeviceName", ""):
        raise ValueError(f"Expected an Apple GPU, found {props.get('DeviceName')!r}")
    if row["elapsed_time"] < job["seconds"]*0.5:
        raise ValueError("Measured interval did not reach benchmark.py's half-target threshold")
    row["ms_per_step"] = 1000*row["elapsed_time"]/row["steps"]
    return row


def write_summary(output, records, tests, repeats, variants=VARIANTS):
    """Keep partial results useful after failures or an interrupted long run."""
    groups = []
    for test in tests:
        for label in variants:
            successful = [record["result"] for record in records if record["test"] == test
                          and record["label"] == label and record["status"] == "passed"]
            groups.append(dict(test=test, label=label, completed=len(successful), expected=repeats,
                median_ns_per_day=statistics.median(row["ns_per_day"] for row in successful) if successful else None,
                median_ms_per_step=statistics.median(row["ms_per_step"] for row in successful) if successful else None))
    save(output / "summary.json", dict(records=records, aggregates=groups))
    lines = ["benchmark.py — single / NVT / HBonds / 4 fs / PME cutoff 0.9 nm",
             "Reported medians are per-case; missing/failed runs are not zero throughput.",
             "test           variant     done       ns/day      ms/step"]
    for group in groups:
        speed = "—" if group["median_ns_per_day"] is None else f"{group['median_ns_per_day']:.3f}"
        step = "—" if group["median_ms_per_step"] is None else f"{group['median_ms_per_step']:.4f}"
        lines.append(f"{group['test']:14} {group['label']:11} {group['completed']}/{repeats:<7} {speed:>10} {step:>12}")
    (output / "summary.txt").write_text("\n".join(lines)+"\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
        allow_abbrev=False,
        epilog="300 seconds is benchmark.py's adaptive target (final interval >=150 seconds), not a timeout.\n"
               "By default: 6 cases x 3 variants x 3 repeats = 54 runs.\n"
               "--variants metal_on selects 18 runs, nominally 1.5 hours plus preparation/calibration.\n"
               "Metal ON includes native float atomics and batching; matrix, hybrid, launch-tuning, and language-target experiments are excluded.\n"
               "Large-force minimization is not exercised by these cases. Equal seeds do not imply\n"
               "identical random streams across platforms. Bindings must be built from this checkout.")
    parser.add_argument("--python", default=sys.executable, help="Python with matching-checkout OpenMM bindings")
    parser.add_argument("--opencl-build", default=str(ROOT / "build/opencl-common-regression"))
    parser.add_argument("--metal-off-build", default=str(ROOT / "build/metal-submission-immediate"))
    parser.add_argument("--metal-on-build", default=str(ROOT / "build/metal-neighbor-optimizations"))
    parser.add_argument("--variants", default=",".join(VARIANTS), help="Comma-separated subset of opencl,metal_off,metal_on")
    parser.add_argument("--tests", default=",".join(TESTS), help="Comma-separated bundled non-plugin cases")
    parser.add_argument("--seconds", type=float, default=300, help="Original benchmark.py target per case (default:300)")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--device", default="0", help="Same single Apple GPU index in all variants")
    parser.add_argument("--output", type=Path, help="New output directory; existing directories are never overwritten")
    parser.add_argument("--jobs", type=int, default=4, help="Build parallelism only; GPU runs are always serial")
    parser.add_argument("--skip-build", action="store_true", help="Use existing binaries; caller must ensure they match cache/source")
    parser.add_argument("--dry-run", action="store_true", help="Read configuration and print plan; no writes/build/OpenMM import")
    parser.add_argument("--check-only", action="store_true", help="Check bindings/library loading, without benchmark integration")
    args = parser.parse_args()
    args.variants = args.variants.split(",")
    if not args.variants or len(set(args.variants)) != len(args.variants) or any(label not in VARIANTS for label in args.variants):
        parser.error(f"Choose distinct variants from {','.join(VARIANTS)}")
    tests = args.tests.split(",")
    if not tests or len(set(tests)) != len(tests) or any(test not in TESTS for test in tests):
        parser.error(f"Choose distinct tests from {','.join(TESTS)}")
    if not math.isfinite(args.seconds) or args.seconds <= 0 or min(args.repeats, args.jobs, args.seed) < 1:
        parser.error("seconds, repeats, jobs, and seed must be positive")
    if args.seed > 2147483647 or args.device != "0":
        parser.error("seed must fit a positive signed int; this single-GPU Metal build requires --device 0")
    for path in (BENCHMARK, *(BENCHMARK.parent / name for name in
                 ("5dfr_minimized.pdb", "5dfr_solv-cube_equil.pdb", "apoa1.pdb"))):
        if not path.is_file():
            parser.error(f"Missing benchmark input: {path}")
    configs = configurations(args)
    labels = [config["label"] for config in configs]
    print(f"Plan: {len(tests)} cases × {len(configs)} variants × {args.repeats} repeats; {args.seconds:g}s target each")
    if "metal_off" in labels:
        print("Metal OFF: fast paths/native float atomics OFF, per-operation submission.")
    if "metal_on" in labels:
        print("Metal ON: fast paths/native float atomics ON, batched submission.")
    for config in configs:
        print(f"  {config['label']}: {config['library']}")
    if args.dry_run:
        return 0
    if sys.platform != "darwin":
        parser.error("This comparison driver requires macOS")
    executable = shutil.which(args.python)
    if executable is None:
        parser.error(f"Python executable not found: {args.python}")
    available = subprocess.run([executable, "-c", "import importlib.util, sys; "
        "sys.exit(0 if importlib.util.find_spec('openmm') else 1)"], capture_output=True, text=True)
    if available.returncode:
        parser.error("Selected Python has no usable OpenMM module. Build/install this checkout's "
                     "Python bindings and select that environment with --python; no build or benchmark was started.")
    output = (args.output or ROOT / "build/benchmarks" / datetime.now().astimezone().strftime("%Y%m%d-%H%M%S")).resolve()
    output.mkdir(parents=True, exist_ok=False)
    (output / "empty-plugins").mkdir()
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    metadata = dict(created=datetime.now(timezone.utc).isoformat(), revision=revision,
        source_status=subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True),
        benchmark_sha256=digest(BENCHMARK), driver_sha256=digest(__file__),
        input_sha256={path.name: digest(path) for path in BENCHMARK.parent.glob("*.pdb")},
        arguments={key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
        configurations=configs, seed_policy="fixed initialization; not cross-platform identical trajectories")
    save(output / "manifest.json", metadata)
    print(f"Output: {output}", flush=True)
    for config in configs:
        if not args.skip_build and not args.check_only:
            code, _ = run_logged(["cmake", "--build", config["build"], "--target",
                                  "OpenMM"+config["platform"], "-j", str(args.jobs)],
                                 output / (config["label"]+"-build.log"))
            if code:
                raise RuntimeError(f"Build failed for {config['label']}; see output directory")
    # A build may have caused CMake to regenerate its cache. Validate and record
    # the resulting settings, not the possibly stale pre-build snapshot.
    configs = configurations(args)
    metadata["configurations"] = configs
    for config in configs:
        config["binary_sha256"] = {key: digest(config[key]) for key in ("core", "library")}
    save(output / "manifest.json", metadata)
    records = []

    def execute(config, test, repeat):
        name = f"{repeat:02d}-{test or 'preflight'}-{config['label']}"
        case = output / name
        case.mkdir()
        job = dict(config=config, test=test, repeat=repeat, seed=args.seed,
                   device=args.device, seconds=args.seconds, revision=revision)
        save(case / "job.json", job)
        environment = dict(os.environ)
        for key in ("MTL_DEBUG_LAYER", "MTL_SHADER_VALIDATION", "MTL_CAPTURE_ENABLED",
                    "DYLD_INSERT_LIBRARIES", "OPENMM_CPU_THREADS"):
            environment.pop(key, None)
        environment.update(OPENMM_PLUGIN_DIR=str(output / "empty-plugins"),
                           DYLD_LIBRARY_PATH=config["build"], PYTHONUNBUFFERED="1")
        record = dict(label=config["label"], test=test, repeat=repeat, exit_code=None,
                      process_wall_seconds=None, directory=str(case), status="running")
        save(case / "status.json", record)
        start = time.monotonic()
        try:
            code, elapsed = run_logged([executable, str(Path(__file__).resolve()), "--worker", str(case / "job.json")],
                                      case / "console.log", environment)
        except BaseException as error:
            record.update(status="interrupted" if isinstance(error, KeyboardInterrupt) else "failed",
                          process_wall_seconds=time.monotonic()-start, error=str(error) or type(error).__name__)
            save(case / "status.json", record)
            if test is not None:
                records.append(record)
            raise
        record.update(exit_code=code, process_wall_seconds=elapsed, status="failed")
        if code:
            record["error"] = f"Worker exited {code}; see console.log"
        else:
            try:
                record["result"] = validate_result(case / "result.json", job) if test else None
                record["status"] = "passed"
            except (ValueError, KeyError, OSError) as error:
                record["error"] = str(error)
        save(case / "status.json", record)
        return record

    # Fail before the long run if bindings or a selected library is unsuitable.
    for config in configs:
        record = execute(config, None, 0)
        if record["status"] != "passed":
            raise RuntimeError(f"Preflight failed: {record['directory']}/console.log")
    if args.check_only:
        print("All library-loading preflights passed; no integration was run.")
        return 0
    write_summary(output, records, tests, args.repeats, labels)
    try:
        for repeat in range(args.repeats):
            offset = repeat % len(configs)
            order = configs[offset:]+configs[:offset]
            for test in tests:
                for config in order:
                    print(f"\n[{repeat+1}/{args.repeats}] {test} / {config['label']}", flush=True)
                    records.append(execute(config, test, repeat+1))
                    write_summary(output, records, tests, args.repeats, labels)
    finally:
        write_summary(output, records, tests, args.repeats, labels)
    print((output / "summary.txt").read_text())
    print(f"Raw data and settings: {output}")
    return int(any(record["status"] != "passed" for record in records))


if __name__ == "__main__":
    try:
        if len(sys.argv) == 3 and sys.argv[1] == "--worker":
            worker(sys.argv[2])
        else:
            sys.exit(main())
    except KeyboardInterrupt:
        print("Interrupted; completed cases and logs are preserved.", file=sys.stderr)
        sys.exit(130)
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
