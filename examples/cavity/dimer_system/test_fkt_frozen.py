#!/usr/bin/env python3
"""
Unit test: F(k,t) must stay constant when positions are frozen.

With the same positions at reference time and at a later time (no motion),
F(k,t) at lag t must equal F(k,0) within numerical tolerance. This guards
against spurious decorrelation from wrong coordinates (e.g. unwrapping drift).
"""

from pathlib import Path
import tempfile
import numpy as np

# Import from same directory
import sys
_here = Path(__file__).resolve().parent
if str(_here) not in sys.path:
    sys.path.insert(0, str(_here))
from fkt_tracker import FKTTracker


def test_fkt_frozen_positions_same_value():
    """F(k,t) with identical positions at t_ref and t_ref+1 ps must equal F(k,0)."""
    np.random.seed(42)
    n_particles = 10
    box_nm = 2.0
    # Positions in [0, box_nm) (wrapped convention)
    positions = np.random.uniform(0, box_nm, size=(n_particles, 3)).astype(np.float64)

    kmag_nm_inv = 113.4
    num_wavevectors = 50
    reference_interval_ps = 100.0  # large so we only get one reference
    max_references = 1
    output_period_ps = 1.0
    t_ref = 1.0

    with tempfile.TemporaryDirectory() as tmpdir:
        prefix = str(Path(tmpdir) / "fkt_frozen_test")
        tracker = FKTTracker(
            kmag_nm_inv=kmag_nm_inv,
            num_wavevectors=num_wavevectors,
            reference_interval_ps=reference_interval_ps,
            max_references=max_references,
            output_period_ps=output_period_ps,
            output_prefix=prefix,
        )
        # Add reference at t_ref
        tracker.update(t_ref, positions)
        # Same positions at t_ref + 1 ps (no motion)
        tracker.update(t_ref + 1.0, positions)
        tracker.finalize()

        # Read output: first line after header is lag=0, second is lag=1
        out_file = Path(prefix + "_fkt_ref_000.txt")
        assert out_file.exists(), f"Expected output file {out_file}"
        lines = out_file.read_text().strip().split("\n")
        data_lines = [l for l in lines if not l.startswith("#") and l.strip()]
        assert len(data_lines) >= 2, f"Expected at least 2 data lines, got {len(data_lines)}"

        # Parse lag_time_ps and F(k,t)
        def parse_row(line):
            parts = line.split()
            return float(parts[0]), float(parts[1])

        lags_and_values = [parse_row(l) for l in data_lines]
        fkt_at_0 = lags_and_values[0][1]
        # Find row with lag ~ 1.0 ps
        row_lag_1 = next((v for lag, v in lags_and_values if abs(lag - 1.0) < 0.01), None)
        assert row_lag_1 is not None, "Expected a row with lag ~ 1.0 ps"
        fkt_at_1 = row_lag_1

    # F(k,t) with same positions must equal F(k,0)
    np.testing.assert_allclose(
        fkt_at_1,
        fkt_at_0,
        rtol=1e-10,
        atol=1e-10,
        err_msg="F(k,t) at lag 1 ps must equal F(k,0) when positions are unchanged (frozen).",
    )


if __name__ == "__main__":
    test_fkt_frozen_positions_same_value()
    print("PASS: test_fkt_frozen_positions_same_value")
