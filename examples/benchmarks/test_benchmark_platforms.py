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

"""Standard-library-only driver tests: no OpenMM imports, builds, or GPU access."""

from contextlib import redirect_stderr, redirect_stdout
import io
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import benchmark_platforms as driver
import run_benchmark_metal_on as wrapper


class BenchmarkVariantsTest(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.args = SimpleNamespace(variants=list(driver.VARIANTS),
            opencl_build=str(self.root / "opencl"), metal_off_build=str(self.root / "off"),
            metal_on_build=str(self.root / "on"))

    def create_cache(self, label):
        platform = "OPENCL" if label == "opencl" else "METAL"
        build = Path(getattr(self.args, label+"_build"))
        build.mkdir()
        cache = dict(CMAKE_BUILD_TYPE="Release", CMAKE_HOME_DIRECTORY=str(driver.ROOT),
                     **{f"OPENMM_BUILD_{platform}_LIB": "ON"})
        if platform == "METAL":
            expected = "ON" if label == "metal_on" else "OFF"
            cache.update(OPENMM_METAL_FAST_BLOCK_BOUNDS=expected,
                         OPENMM_METAL_NATIVE_FLOAT_ATOMICS=expected,
                         OPENMM_METAL_RECORD_AND_COMMIT="OFF" if label == "metal_on" else "ON",
                         OPENMM_METAL_EXPERIMENTAL_MATRIX_SCREEN="OFF")
        (build / "CMakeCache.txt").write_text("\n".join(f"{key}:STRING={value}" for key, value in cache.items()))
        return build

    def test_only_selected_cache_is_required(self):
        self.create_cache("metal_on")
        self.args.variants = ["metal_on"]
        configs = driver.configurations(self.args)
        self.assertEqual([config["label"] for config in configs], ["metal_on"])

    def test_default_includes_all_three_variants(self):
        for label in driver.VARIANTS:
            self.create_cache(label)
        self.assertEqual([config["label"] for config in driver.configurations(self.args)], list(driver.VARIANTS))

    def test_both_metal_variants_still_require_same_switches(self):
        for label in driver.VARIANTS:
            self.create_cache(label)
        cache = Path(self.args.metal_on_build) / "CMakeCache.txt"
        cache.write_text(cache.read_text()+"\nOPENMM_METAL_FAST_ADDITIONAL:BOOL=ON\n")
        with self.assertRaisesRegex(ValueError, "different fast-path switches"):
            driver.configurations(self.args)

    def test_on_only_still_rejects_disabled_optimization(self):
        build = self.create_cache("metal_on")
        self.args.variants = ["metal_on"]
        cache = build / "CMakeCache.txt"
        cache.write_text(cache.read_text().replace("OPENMM_METAL_FAST_BLOCK_BOUNDS:STRING=ON",
                                                 "OPENMM_METAL_FAST_BLOCK_BOUNDS:STRING=OFF"))
        with self.assertRaisesRegex(ValueError, "wrong settings"):
            driver.configurations(self.args)

    def test_hybrid_experiment_cannot_masquerade_as_standard_endpoint(self):
        for label in ("metal_off", "metal_on"):
            with self.subTest(label=label):
                build = self.create_cache(label)
                self.args.variants = [label]
                cache = build / "CMakeCache.txt"
                original = cache.read_text()
                cache.write_text(original+"\nOPENMM_METAL_EXPERIMENTAL_NONBONDED_HYBRID:BOOL=OFF\n")
                self.assertEqual(driver.configurations(self.args)[0]["label"], label)
                cache.write_text(original+"\nOPENMM_METAL_EXPERIMENTAL_NONBONDED_HYBRID:BOOL=ON\n")
                with self.assertRaisesRegex(ValueError, "EXPERIMENTAL_NONBONDED_HYBRID"):
                    driver.configurations(self.args)

    def test_tuning_cannot_masquerade_as_standard_endpoint(self):
        switches = ("OPENMM_METAL_TUNE_FORCE_THREADGROUP_SIZE",
                    "OPENMM_METAL_TUNE_FORCE_GROUPS_PER_COMPUTE_UNIT",
                    "OPENMM_METAL_TUNE_FORCE_PIPELINE_MAX_THREADS",
                    "OPENMM_METAL_TUNE_FORCE_REQUIRED_THREADS",
                    "OPENMM_METAL_TUNE_FORCE_FUTURE_OPTION",
                    "OPENMM_METAL_TUNE_LANGUAGE_VERSION")
        for label in ("metal_off", "metal_on"):
            build = self.create_cache(label)
            self.args.variants = [label]
            cache = build / "CMakeCache.txt"
            original = cache.read_text()
            for switch in switches:
                with self.subTest(label=label, switch=switch):
                    cache.write_text(original+f"\n{switch}:BOOL=OFF\n")
                    self.assertEqual(driver.configurations(self.args)[0]["label"], label)
                    cache.write_text(original+f"\n{switch}:BOOL=ON\n")
                    with self.assertRaisesRegex(ValueError, switch):
                        driver.configurations(self.args)

    def test_disabled_launch_tuning_preserves_numeric_values_in_cmake_metadata(self):
        settings = dict(OPENMM_METAL_FORCE_THREADGROUP_SIZE="64",
                        OPENMM_METAL_FORCE_GROUPS_PER_COMPUTE_UNIT="3",
                        OPENMM_METAL_FORCE_PIPELINE_MAX_THREADS="128")
        for label in ("metal_off", "metal_on"):
            with self.subTest(label=label):
                build = self.create_cache(label)
                self.args.variants = [label]
                cache = build / "CMakeCache.txt"
                cache.write_text(cache.read_text()+"\n"+"\n".join(
                    f"{key}:STRING={value}\n{key.replace('OPENMM_METAL_FORCE_', 'OPENMM_METAL_TUNE_FORCE_')}:BOOL=OFF"
                    for key, value in settings.items()))
                config = driver.configurations(self.args)[0]
                metadata = self.root / (label+"-metadata.json")
                driver.save(metadata, dict(configurations=[config]))
                saved = json.loads(metadata.read_text())["configurations"][0]
                for key, value in settings.items():
                    self.assertEqual(saved["cmake"][key], value)
                    self.assertNotIn(key, saved["flags"])
                    self.assertEqual(saved["flags"][key.replace("OPENMM_METAL_FORCE_", "OPENMM_METAL_TUNE_FORCE_")], "OFF")

    def test_disabled_language_override_preserves_target_in_metadata(self):
        for label in ("metal_off", "metal_on"):
            with self.subTest(label=label):
                build = self.create_cache(label)
                self.args.variants = [label]
                cache = build / "CMakeCache.txt"
                cache.write_text(cache.read_text()+"\nOPENMM_METAL_TUNE_LANGUAGE_VERSION:BOOL=OFF"
                                 "\nOPENMM_METAL_LANGUAGE_VERSION:STRING=410\n")
                config = driver.configurations(self.args)[0]
                self.assertEqual(config["cmake"]["OPENMM_METAL_LANGUAGE_VERSION"], "410")
                self.assertNotIn("OPENMM_METAL_LANGUAGE_VERSION", config["flags"])
                self.assertEqual(config["flags"]["OPENMM_METAL_TUNE_LANGUAGE_VERSION"], "OFF")

    def test_summary_contains_only_selected_variant(self):
        records = [dict(test="pme", label="metal_on", status="passed",
                        result=dict(ns_per_day=speed, ms_per_step=1000/speed)) for speed in (10, 20, 30)]
        driver.write_summary(self.root, records, ["pme"], 3, ["metal_on"])
        summary = json.loads((self.root / "summary.json").read_text())
        self.assertEqual(len(summary["aggregates"]), 1)
        self.assertEqual(summary["aggregates"][0]["median_ns_per_day"], 20)
        self.assertEqual(summary["aggregates"][0]["completed"], 3)
        self.assertNotIn("metal_off", (self.root / "summary.txt").read_text())
        self.assertNotIn("opencl", (self.root / "summary.txt").read_text())

    def test_default_summary_keeps_all_variants(self):
        driver.write_summary(self.root, [], driver.TESTS, 3)
        summary = json.loads((self.root / "summary.json").read_text())
        self.assertEqual(len(summary["aggregates"]), 18)

    def test_dry_run_does_not_build_import_openmm_or_create_output(self):
        build = self.create_cache("metal_on")
        output = self.root / "results"
        argv = [str(driver.__file__), "--variants", "metal_on", "--metal-on-build", str(build),
                "--opencl-build", str(self.root / "missing-opencl"),
                "--metal-off-build", str(self.root / "missing-off"), "--output", str(output), "--dry-run"]
        console = io.StringIO()
        with patch("sys.argv", argv), patch.object(driver, "run_logged") as run, redirect_stdout(console):
            self.assertEqual(driver.main(), 0)
        run.assert_not_called()
        self.assertFalse(output.exists())
        self.assertIn("6 cases × 1 variants × 3 repeats; 300s target each", console.getvalue())

    def test_invalid_variants_are_rejected(self):
        for variants in ("", "metal_on,metal_on", "unknown"):
            with self.subTest(variants=variants), patch("sys.argv", [str(driver.__file__), "--variants", variants]), redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as error:
                    driver.main()
                self.assertEqual(error.exception.code, 2)

    def test_wrapper_uses_prepared_python_and_forwards_options(self):
        with patch("sys.argv", [str(wrapper.__file__), "--dry-run"]), patch.object(Path, "is_file", return_value=True), patch.object(wrapper.os, "execv") as execute:
            wrapper.main()
        python = str(wrapper.ROOT / "build/benchmark-venv/bin/python")
        execute.assert_called_once_with(python, [python, str(wrapper.ROOT / "examples/benchmarks/benchmark_platforms.py"),
                                               "--variants", "metal_on", "--skip-build", "--dry-run"])

    def test_wrapper_prevents_variant_override(self):
        for option in (["--variants", "opencl"], ["--variants=opencl"]):
            with self.subTest(option=option), patch("sys.argv", [str(wrapper.__file__), *option]):
                with self.assertRaisesRegex(SystemExit, "only metal_on"):
                    wrapper.main()


if __name__ == "__main__":
    unittest.main()
