#!/usr/bin/env python3
"""Poll i-PI / UMA driver logs every 60 s while a run is in progress.

Usage (from tests/uma_ice_rpmd/, in a second terminal while orchestrator runs):

  python monitor_ipi_uma_run.py

Optional: ``--interval 60`` (seconds), ``--minutes 120`` max wall time.

Exits non-zero if the i-PI log shows a traceback or known fatal strings.
Does not start or stop the simulation — read-only monitoring.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

_SCRIPT_DIR = Path(__file__).resolve().parent
IPI_LOG = _SCRIPT_DIR / "ipi" / "i-pi_run.log"
DRIVER_LOG = _SCRIPT_DIR / "pipeline_out" / "uma_driver.log"

_FATAL = (
    "Traceback (most recent call last)",
    "OutOfMemoryError",
    "CUDA out of memory",
    "Address already in use",
    "ValueError: cannot reshape array of size 0",
)


def _tail(path: Path, n: int = 25) -> str:
    if not path.is_file():
        return f"(no file: {path})\n"
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    return "\n".join(lines[-n:]) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser(description="Monitor i-PI + UMA driver logs")
    ap.add_argument("--interval", type=float, default=60.0, help="Seconds between checks")
    ap.add_argument("--minutes", type=float, default=0.0, help="Stop after this many minutes (0 = unlimited)")
    ap.add_argument("--ipi-log", type=Path, default=IPI_LOG)
    ap.add_argument("--driver-log", type=Path, default=DRIVER_LOG)
    args = ap.parse_args()

    t0 = time.time()
    n = 0
    while True:
        n += 1
        elapsed = int(time.time() - t0)
        print(f"\n{'=' * 60}\n# check {n}  (+{elapsed}s)\n{'=' * 60}", flush=True)

        # lightweight process hint
        try:
            out = subprocess.run(
                ["pgrep", "-af", "i-pi|i-pi-py_driver|run_ipi_lammps_uma_rpmd"],
                capture_output=True,
                text=True,
                timeout=5,
            )
            procs = out.stdout.strip() or "(no matching processes)"
            print("Processes:\n", procs[:2000], "\n", flush=True)
        except (OSError, subprocess.TimeoutExpired) as e:
            print(f"(pgrep unavailable: {e})", flush=True)

        ipi_txt = args.ipi_log.read_text(encoding="utf-8", errors="replace") if args.ipi_log.is_file() else ""
        for sig in _FATAL:
            if sig in ipi_txt:
                print(f"FATAL in i-PI log: matched {sig!r}", file=sys.stderr, flush=True)
                print(_tail(args.ipi_log, 40), flush=True)
                sys.exit(1)

        print("--- i-pi_run.log (last 20 lines) ---", flush=True)
        print(_tail(args.ipi_log, 20), flush=True)
        print("--- uma_driver.log (last 15 lines) ---", flush=True)
        print(_tail(args.driver_log, 15), flush=True)

        if args.minutes > 0 and (time.time() - t0) >= args.minutes * 60:
            print("Max minutes reached; exiting monitor.", flush=True)
            break
        time.sleep(args.interval)


if __name__ == "__main__":
    main()
