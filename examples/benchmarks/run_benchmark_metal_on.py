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

"""Run only the existing Metal ON build with the prepared benchmark environment.

Defaults: six benchmark.py cases, three repeats, 300-second adaptive target.
Results go to a new timestamped build/benchmarks directory. No build or package
installation is performed. Additional benchmark_platforms.py options such as
--dry-run, --check-only, --output, --seconds, and --repeats are forwarded.
"""

import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]


def main():
    """Use checkout-relative paths regardless of the caller's working directory."""
    if any(argument == "--variants" or argument.startswith("--variants=") for argument in sys.argv[1:]):
        raise SystemExit("This entry point runs only metal_on; use benchmark_platforms.py to select other variants.")
    python = ROOT / "build/benchmark-venv/bin/python"
    if not python.is_file():
        raise SystemExit(f"Prepared benchmark Python was not found: {python}")
    driver = ROOT / "examples/benchmarks/benchmark_platforms.py"
    os.execv(str(python), [str(python), str(driver), "--variants", "metal_on", "--skip-build", *sys.argv[1:]])


if __name__ == "__main__":
    main()
