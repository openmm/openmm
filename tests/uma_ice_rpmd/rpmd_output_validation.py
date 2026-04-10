#!/usr/bin/env python3
"""Validate RPMD order CSV completeness for pipeline runs."""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class OrderCsvSummary:
    """Minimal summary of an RPMD order CSV."""

    row_count: int
    first_time_ps: float | None
    last_time_ps: float | None


def summarize_order_csv(path: Path) -> OrderCsvSummary:
    """Return row count and first/last ``time_ps`` from an order CSV."""
    row_count = 0
    first_time_ps: float | None = None
    last_time_ps: float | None = None

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or "time_ps" not in reader.fieldnames:
            raise ValueError(f"{path} is missing required 'time_ps' column")
        for row in reader:
            time_ps = float(row["time_ps"])
            if first_time_ps is None:
                first_time_ps = time_ps
            last_time_ps = time_ps
            row_count += 1

    return OrderCsvSummary(
        row_count=row_count,
        first_time_ps=first_time_ps,
        last_time_ps=last_time_ps,
    )


def validate_order_csv_complete(
    path: Path,
    *,
    expected_final_time_ps: float,
    label: str = "order CSV",
    tolerance_ps: float = 1.0e-9,
) -> OrderCsvSummary:
    """Validate that an order CSV exists and reaches the requested final time."""
    if not path.is_file():
        raise ValueError(f"{label} missing: {path}")

    summary = summarize_order_csv(path)
    if summary.row_count == 0 or summary.last_time_ps is None:
        raise ValueError(f"{label} has no data rows: {path}")
    if summary.last_time_ps + tolerance_ps < expected_final_time_ps:
        raise ValueError(
            f"{label} incomplete: expected final time >= {expected_final_time_ps:.6f} ps, "
            f"found {summary.last_time_ps:.6f} ps in {path}"
        )
    return summary


def main() -> None:
    ap = argparse.ArgumentParser(description="Validate RPMD order CSV completeness")
    ap.add_argument("--csv", type=Path, required=True, help="CSV to validate")
    ap.add_argument(
        "--expected-final-ps",
        type=float,
        required=True,
        help="Required final time_ps for the CSV to be considered complete",
    )
    ap.add_argument("--label", type=str, default="order CSV")
    args = ap.parse_args()

    summary = validate_order_csv_complete(
        args.csv,
        expected_final_time_ps=args.expected_final_ps,
        label=args.label,
    )
    print(
        f"{args.label} complete: rows={summary.row_count} "
        f"time_ps={summary.first_time_ps:.6f}->{summary.last_time_ps:.6f}"
    )


if __name__ == "__main__":
    main()
