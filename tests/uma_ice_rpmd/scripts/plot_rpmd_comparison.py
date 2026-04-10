#!/usr/bin/env python3
"""
Plot ice **order-parameter CVs** (⟨Q6⟩, ⟨q_tet⟩) and kinetic T from RPMD exports.

The structural collective variables (CVs) tracked here are **Steinhardt Q6** and
**Chau–Harding q_tet**, computed with the path-integral / bead-averaging estimator
(see ``ice_order_parameters.ice_order_metrics_path_integral``). The default figure
shows **time series** of these CVs plus centroid kinetic temperature.

Use ``--cv-plane PATH`` to also save the **trajectory projected into the
(⟨Q6⟩, ⟨q_tet⟩) plane** (same data as the top two panels, useful for visualizing
sampling in CV space).

Third panel: ``T_K`` (centroid kinetic temperature) vs time, with optional bath reference
line (default 243 K for PILE-G). Use this to check thermostat behaviour.

Reads CSV files produced by:
  - OpenMM UMA RPMD: test_uma_ice_rpmd.py (ice_order_openmm_rpmd.csv)
  - TIP4P RPMD: run_openmm_tip4p_rpmd.py (ice_order_tip4p_rpmd.csv)
  - i-PI+LAMMPS: ipi_order_from_traj.py (ice_order_ipi_rpmd.csv; default per-atom PBC wrap)

Usage:
  cd tests/uma_ice_rpmd
  python plot_rpmd_comparison.py --openmm pipeline_out/ice_order_openmm_rpmd.csv --ipi pipeline_out/ice_order_ipi_rpmd.csv --tip4p pipeline_out/ice_order_tip4p_rpmd.csv -o rpmd_comparison.png

  # Same inputs, plus CV-plane (Q6 vs q_tet) trajectory:
  python plot_rpmd_comparison.py ... -o rpmd_comparison.png --cv-plane pipeline_out/rpmd_cv_plane.png

  # Three RPMD curves (unstable CUDA vs stable Reference CPU integrator vs i-PI+LAMMPS):
  python plot_rpmd_comparison.py --openmm pipeline_out/ice_order_openmm_rpmd_1ps.csv \\
    --openmm-reference pipeline_out/reference_platform_test.csv \\
    --ipi pipeline_out/ice_order_ipi_rpmd.csv -o pipeline_out/rpmd_comparison.png \\
    --cv-plane pipeline_out/rpmd_cv_plane.png

  # Time window (default ``--time-max-ps`` is 0.1); full span: ``--time-max-ps 0``.

For **classical (single-bead) MD** only — LAMMPS vs OpenMM vs TIP4P — use ``plot_md_comparison.py`` (not this script).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _HAS_MPL = True
except ImportError:
    _HAS_MPL = False

_ROOT = Path(__file__).resolve().parent.parent

# Compact figure + larger type; serif + Computer Modern mathtext ≈ LaTeX default look
# (no external ``latex`` required; see ``mathtext.fontset``).
_PLOT_RC = {
    "font.size": 12,
    "font.family": "serif",
    "font.serif": ["DejaVu Serif", "Bitstream Vera Serif", "Computer Modern Roman", "Times", "Palatino"],
    "mathtext.fontset": "cm",
    "axes.labelsize": 14,
    "axes.titlesize": 13,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
    "figure.titlesize": 14,
    "axes.unicode_minus": False,
}
# Width fixed; height gives room for 3 stacked panels without a huge title gap.
_FIGSIZE_IN = (6.0, 7.75)
_DEFAULT_DPI = 150


def _format_title_two_lines(title: str) -> str:
    """Return a two-line suptitle when possible (shorter single line per row).

    If *title* already contains a newline, it is returned unchanged. Otherwise,
    split on an em dash `` — `` (common in our CLI titles) into headline + subtitle.
    """
    t = title.strip()
    if "\n" in t:
        return t
    em = " — "
    if em in t:
        left, right = t.split(em, 1)
        return f"{left.rstrip()}\n{right.strip()}"
    return t


def _read_n_oxygen_from_csv(path: Path) -> int | None:
    """Return ``n_oxygen`` from the first parsable data row, or ``None`` if absent.

    For water ice order exports, ``n_oxygen`` equals the number of water molecules.
    """
    if not path.is_file():
        return None
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 2:
        return None
    cols = [c.strip() for c in lines[0].split(",")]
    if "n_oxygen" not in cols:
        return None
    idx = cols.index("n_oxygen")
    for line in lines[1:]:
        parts = [p.strip() for p in line.split(",")]
        if len(parts) <= idx:
            continue
        try:
            return int(float(parts[idx]))
        except ValueError:
            continue
    return None


def _enrich_title_with_system(
    title: str,
    n_molecules: int | None,
    n_beads: int | None,
) -> str:
    """Append `` · N molecules, P beads`` to the last title line when values are set."""
    parts: list[str] = []
    if n_molecules is not None:
        parts.append(f"{n_molecules} molecules")
    if n_beads is not None:
        parts.append(f"{n_beads} beads")
    if not parts:
        return title
    suffix = ", ".join(parts)
    lines = title.split("\n")
    lines[-1] = f"{lines[-1]} · {suffix}"
    return "\n".join(lines)


def _read_csv(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray | None]:
    """Return (time_ps, q6_mean, q_tet_mean, T_K) from ice order CSV.

    ``T_K`` is **centroid kinetic temperature** (same estimator as OpenMM
    ``centroid_kinetic_energy_and_temperature`` / i-PI post-process from bead
    velocities when vtraj is available). Missing or empty cells become ``nan``. If the
    ``T_K`` column is absent, the fourth return is ``None`` (e.g. some i-PI exports).
    """
    lines = path.read_text().splitlines()
    if len(lines) < 2:
        return np.array([]), np.array([]), np.array([]), None
    header = lines[0]
    cols = [c.strip() for c in header.split(",")]
    idx_t = cols.index("time_ps")
    idx_q6 = cols.index("q6_mean")
    idx_qt = cols.index("q_tet_mean")
    has_tk = "T_K" in cols
    idx_tk = cols.index("T_K") if has_tk else -1
    t_list, q6_list, qt_list = [], [], []
    tk_list: list[float] = []
    for line in lines[1:]:
        parts = [p.strip() for p in line.split(",")]
        if len(parts) <= max(idx_t, idx_q6, idx_qt):
            continue
        try:
            t_list.append(float(parts[idx_t]))
            q6_list.append(float(parts[idx_q6]))
            qt_list.append(float(parts[idx_qt]))
        except ValueError:
            continue
        if has_tk:
            if idx_tk < len(parts) and parts[idx_tk] != "":
                try:
                    tk_list.append(float(parts[idx_tk]))
                except ValueError:
                    tk_list.append(float("nan"))
            else:
                tk_list.append(float("nan"))
    t_arr = np.array(t_list)
    q6_arr = np.array(q6_list)
    qt_arr = np.array(qt_list)
    tk_arr: np.ndarray | None = np.array(tk_list, dtype=np.float64) if has_tk else None
    if tk_arr is not None and tk_arr.shape != t_arr.shape:
        # Should not happen if every parsed row appends tk
        tk_arr = None
    return t_arr, q6_arr, qt_arr, tk_arr


def _apply_time_max_ps(
    t_max_ps: float | None,
    t_arr: np.ndarray,
    q6: np.ndarray,
    qt: np.ndarray,
    tk: np.ndarray | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray | None]:
    """Return arrays restricted to ``0 <= time_ps <= t_max_ps`` when *t_max_ps* > 0."""
    if t_max_ps is None or t_max_ps <= 0 or t_arr.size == 0:
        return t_arr, q6, qt, tk
    mask = (t_arr >= 0.0) & (t_arr <= float(t_max_ps))
    if not np.any(mask):
        return t_arr, q6, qt, tk
    tk_out = tk[mask] if tk is not None else None
    return t_arr[mask], q6[mask], qt[mask], tk_out


def _save_cv_plane_figure(
    output: Path,
    *,
    has_openmm: bool,
    q6_o: np.ndarray,
    qt_o: np.ndarray,
    reference_data: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray | None] | None,
    cuda_sequential_data: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray | None] | None,
    has_ipi: bool,
    q6_i: np.ndarray,
    qt_i: np.ndarray,
    tip4p_data: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray | None] | None,
    title: str,
    n_molecules: int | None,
    n_beads_ann: int | None,
    dpi: int,
) -> None:
    """Save ⟨Q6⟩ vs ⟨q_tet⟩ (RPMD order-parameter CV plane); points connected in time order."""
    with plt.rc_context(_PLOT_RC):
        fig, ax = plt.subplots(1, 1, figsize=(6.0, 5.0), constrained_layout=True)
        if has_openmm and q6_o.size:
            ax.plot(q6_o, qt_o, "b-", label="OpenMM CUDA + UMA", alpha=0.8, lw=1.8)
            ax.scatter(q6_o[0], qt_o[0], c="blue", s=36, marker="o", zorder=5, edgecolors="k", linewidths=0.4)
            ax.scatter(q6_o[-1], qt_o[-1], c="blue", s=36, marker="s", zorder=5, edgecolors="k", linewidths=0.4)
        if reference_data is not None:
            _, q6_r, qt_r, _ = reference_data
            if q6_r.size:
                ax.plot(q6_r, qt_r, "c-", label="OpenMM Reference + UMA", alpha=0.85, lw=1.8)
                ax.scatter(q6_r[0], qt_r[0], c="cyan", s=36, marker="o", zorder=5, edgecolors="k", linewidths=0.4)
                ax.scatter(q6_r[-1], qt_r[-1], c="cyan", s=36, marker="s", zorder=5, edgecolors="k", linewidths=0.4)
        if cuda_sequential_data is not None:
            _, q6_s, qt_s, _ = cuda_sequential_data
            if q6_s.size:
                ax.plot(
                    q6_s,
                    qt_s,
                    "m-.",
                    label="OpenMM CUDA + UMA (sequential)",
                    alpha=0.85,
                    lw=1.6,
                )
                ax.scatter(q6_s[0], qt_s[0], c="magenta", s=36, marker="o", zorder=5, edgecolors="k", linewidths=0.4)
                ax.scatter(q6_s[-1], qt_s[-1], c="magenta", s=36, marker="s", zorder=5, edgecolors="k", linewidths=0.4)
        if has_ipi and q6_i.size:
            ax.plot(q6_i, qt_i, "r--", label="i-PI + LAMMPS UMA", alpha=0.8, lw=1.8)
            ax.scatter(q6_i[0], qt_i[0], c="red", s=36, marker="o", zorder=5, edgecolors="k", linewidths=0.4)
            ax.scatter(q6_i[-1], qt_i[-1], c="red", s=36, marker="s", zorder=5, edgecolors="k", linewidths=0.4)
        if tip4p_data is not None:
            _, q6_t, qt_t, _ = tip4p_data
            if q6_t.size:
                ax.plot(q6_t, qt_t, "g-.", label="TIP4P/2005f RPMD", alpha=0.85, lw=1.8)
                ax.scatter(q6_t[0], qt_t[0], c="green", s=36, marker="o", zorder=5, edgecolors="k", linewidths=0.4)
                ax.scatter(q6_t[-1], qt_t[-1], c="green", s=36, marker="s", zorder=5, edgecolors="k", linewidths=0.4)
        ax.set_xlabel(r"$\langle Q_6 \rangle$ (order CV)")
        ax.set_ylabel(r"$\langle q_{\mathrm{tet}} \rangle$ (order CV)")
        ax.set_title("RPMD trajectory in order-parameter plane")
        ax.legend()
        ax.grid(True, alpha=0.3)
        cv_title = _format_title_two_lines(title) + "\n(circle=start, square=end)"
        fig.suptitle(_enrich_title_with_system(cv_title, n_molecules, n_beads_ann))
        output.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output, dpi=dpi, bbox_inches="tight", pad_inches=0.08)


def main() -> None:
    ap = argparse.ArgumentParser(description="Plot OpenMM vs i-PI+LAMMPS Q6 comparison")
    ap.add_argument(
        "--openmm",
        type=Path,
        default=_ROOT / "pipeline_out" / "ice_order_openmm_rpmd.csv",
        help="OpenMM UMA RPMD order CSV (test_uma_ice_rpmd / run_openmm_rpmd_reference).",
    )
    ap.add_argument(
        "--openmm-reference",
        type=Path,
        default=None,
        metavar="PATH",
        help=(
            "Optional second OpenMM curve (same CSV schema), e.g. Reference CPU integrator "
            "with UMA forces — plotted cyan for comparison to --openmm (often CUDA)."
        ),
    )
    ap.add_argument(
        "--openmm-sequential",
        type=Path,
        default=None,
        metavar="PATH",
        help=(
            "Optional third OpenMM CUDA curve: same schema as --openmm, e.g. "
            "OPENMMML_RPMD_NO_BATCH=1 sequential per-bead forces — plotted magenta."
        ),
    )
    ap.add_argument(
        "--ipi",
        type=Path,
        default=None,
        help=(
            "i-PI+LAMMPS order CSV; omit to skip i-PI curves "
            f"(default pipeline path: {_ROOT / 'pipeline_out' / 'ice_order_ipi_rpmd.csv'})."
        ),
    )
    ap.add_argument(
        "--tip4p",
        type=Path,
        default=None,
        help="TIP4P/2005f OpenMM RPMD order CSV (run_openmm_tip4p_rpmd.py)",
    )
    ap.add_argument("-o", "--output", type=Path, default=_ROOT / "pipeline_out" / "rpmd_comparison.png")
    ap.add_argument(
        "--cv-plane",
        type=Path,
        default=None,
        metavar="PATH",
        help=(
            "Also save a figure: ⟨Q6⟩ vs ⟨q_tet⟩ (RPMD order-parameter CV plane), "
            "trajectory in time order (circle=start, square=end)."
        ),
    )
    ap.add_argument(
        "--title",
        type=str,
        default=(
            "RPMD ice order: OpenMM UMA vs i-PI+LAMMPS UMA\n"
            "(243 K, PILE-G)"
        ),
        help="Figure suptitle; use \\n for manual line breaks, or rely on splitting at ' — '.",
    )
    ap.add_argument(
        "--bath-temperature-k",
        type=float,
        default=243.0,
        help="Target bath temperature (K) for horizontal reference line on the T_K panel (PILE-G).",
    )
    ap.add_argument(
        "--tk-ymax-factor",
        type=float,
        default=1.1,
        metavar="F",
        help=(
            "T_K panel only: y-axis top = F × bath temperature (default 1.1 = bath +10%%). "
            "Data above this are clipped from view."
        ),
    )
    ap.add_argument(
        "--no-temperature-line",
        action="store_true",
        help="Do not draw the bath-temperature reference line on the kinetic T panel.",
    )
    ap.add_argument(
        "--dpi",
        type=int,
        default=_DEFAULT_DPI,
        help=f"Raster resolution for PNG output (default {_DEFAULT_DPI}).",
    )
    ap.add_argument(
        "--molecules",
        type=int,
        default=None,
        metavar="N",
        help="Molecule count in suptitle; default: read n_oxygen from the first available CSV.",
    )
    ap.add_argument(
        "--beads",
        type=int,
        default=32,
        metavar="P",
        help=(
            "RPMD bead count in suptitle (default 32, matching run_openmm_uma_rpmd_only.sh). "
            "Use 0 to omit beads from the annotation."
        ),
    )
    ap.add_argument(
        "--time-max-ps",
        type=float,
        default=0.1,
        metavar="T",
        help=(
            "Clip all time-series curves to 0…T ps and fix the x-axis to [0, T] (default 0.1). "
            "Use 0 or a negative value to plot each series over its full recorded time range."
        ),
    )
    args = ap.parse_args()

    if not _HAS_MPL:
        print("matplotlib required: pip install matplotlib", file=sys.stderr)
        sys.exit(1)

    has_openmm = args.openmm.is_file()
    ipi_path = args.ipi
    has_ipi = ipi_path is not None and Path(ipi_path).is_file()
    tip4p_data: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray | None] | None = None
    if args.tip4p is not None and args.tip4p.is_file():
        tip4p_data = _read_csv(args.tip4p)

    ref_path = args.openmm_reference
    reference_data: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray | None] | None = None
    if ref_path is not None and Path(ref_path).is_file():
        reference_data = _read_csv(Path(ref_path))

    seq_path = args.openmm_sequential
    cuda_sequential_data: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray | None] | None = None
    if seq_path is not None and Path(seq_path).is_file():
        cuda_sequential_data = _read_csv(Path(seq_path))

    if (
        not has_openmm
        and not has_ipi
        and tip4p_data is None
        and reference_data is None
        and cuda_sequential_data is None
    ):
        print(
            "Need at least one of: existing --openmm CSV, existing --ipi CSV, "
            "--openmm-reference CSV, --openmm-sequential CSV, or --tip4p CSV.",
            file=sys.stderr,
        )
        sys.exit(1)

    t_o = q6_o = qt_o = np.array([])
    tk_o: np.ndarray | None = None
    if has_openmm:
        t_o, q6_o, qt_o, tk_o = _read_csv(args.openmm)
    t_i = q6_i = qt_i = np.array([])
    tk_i: np.ndarray | None = None
    if has_ipi:
        t_i, q6_i, qt_i, tk_i = _read_csv(Path(ipi_path))

    t_max = args.time_max_ps if args.time_max_ps > 0 else None
    if has_openmm:
        t_o, q6_o, qt_o, tk_o = _apply_time_max_ps(t_max, t_o, q6_o, qt_o, tk_o)
    if has_ipi:
        t_i, q6_i, qt_i, tk_i = _apply_time_max_ps(t_max, t_i, q6_i, qt_i, tk_i)
    if reference_data is not None:
        t_r, q6_r, qt_r, tk_r = reference_data
        t_r, q6_r, qt_r, tk_r = _apply_time_max_ps(t_max, t_r, q6_r, qt_r, tk_r)
        reference_data = (t_r, q6_r, qt_r, tk_r)
    if cuda_sequential_data is not None:
        t_s, q6_s, qt_s, tk_s = cuda_sequential_data
        t_s, q6_s, qt_s, tk_s = _apply_time_max_ps(t_max, t_s, q6_s, qt_s, tk_s)
        cuda_sequential_data = (t_s, q6_s, qt_s, tk_s)
    if tip4p_data is not None:
        t_t, q6_t, qt_t, tk_t = tip4p_data
        t_t, q6_t, qt_t, tk_t = _apply_time_max_ps(t_max, t_t, q6_t, qt_t, tk_t)
        tip4p_data = (t_t, q6_t, qt_t, tk_t)

    def _plot_tk(
        ax: "plt.Axes",
        t_arr: np.ndarray,
        tk_arr: np.ndarray | None,
        label: str,
        style: str,
        lw: float | None = None,
    ) -> None:
        if tk_arr is None or t_arr.size == 0:
            return
        mask = np.isfinite(tk_arr)
        if not np.any(mask):
            return
        ax.plot(t_arr[mask], tk_arr[mask], style, label=label, alpha=0.85, lw=lw)

    title_fmt = _format_title_two_lines(args.title)
    n_mol = args.molecules
    if n_mol is None:
        if has_openmm:
            n_mol = _read_n_oxygen_from_csv(args.openmm)
        if n_mol is None and args.tip4p is not None:
            n_mol = _read_n_oxygen_from_csv(args.tip4p)
        if n_mol is None and has_ipi and ipi_path is not None:
            n_mol = _read_n_oxygen_from_csv(Path(ipi_path))
        if n_mol is None and ref_path is not None:
            n_mol = _read_n_oxygen_from_csv(Path(ref_path))
        if n_mol is None and seq_path is not None:
            n_mol = _read_n_oxygen_from_csv(Path(seq_path))
    n_beads_ann: int | None = args.beads if args.beads > 0 else None

    with plt.rc_context(_PLOT_RC):
        # constrained_layout pulls panels up under a multi-line suptitle (less dead space
        # than tight_layout(rect=...) + bbox_inches=tight, which often over-reserves the top).
        fig, (ax1, ax2, ax3) = plt.subplots(
            3,
            1,
            figsize=_FIGSIZE_IN,
            sharex=True,
            constrained_layout=True,
        )
        if has_openmm and t_o.size:
            ax1.plot(t_o, q6_o, "b-", label="OpenMM CUDA + UMA", alpha=0.8)
            ax2.plot(t_o, qt_o, "b-", label="OpenMM CUDA + UMA", alpha=0.8)
        if reference_data is not None:
            t_r, q6_r, qt_r, tk_r = reference_data
            if t_r.size:
                ax1.plot(t_r, q6_r, "c-", label="OpenMM Reference + UMA", alpha=0.85)
                ax2.plot(t_r, qt_r, "c-", label="OpenMM Reference + UMA", alpha=0.85)
        if cuda_sequential_data is not None:
            t_s, q6_s, qt_s, _ = cuda_sequential_data
            if t_s.size:
                ax1.plot(
                    t_s,
                    q6_s,
                    "m-.",
                    label="OpenMM CUDA + UMA (sequential)",
                    alpha=0.85,
                    lw=1.6,
                )
                ax2.plot(t_s, qt_s, "m-.", label="OpenMM CUDA + UMA (sequential)", alpha=0.85, lw=1.6)
        if has_ipi and t_i.size:
            ax1.plot(t_i, q6_i, "r--", label="i-PI + LAMMPS UMA", alpha=0.8)
            ax2.plot(t_i, qt_i, "r--", label="i-PI + LAMMPS UMA", alpha=0.8)
        if tip4p_data is not None:
            t_t, q6_t, qt_t, tk_t = tip4p_data
            ax1.plot(t_t, q6_t, "g-.", label="TIP4P/2005f RPMD", alpha=0.85, lw=1.8)
            ax2.plot(t_t, qt_t, "g-.", label="TIP4P/2005f RPMD", alpha=0.85, lw=1.8)
        ax1.set_ylabel(r"$\langle Q_6 \rangle$")
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        ax2.set_ylabel(r"$\langle q_{\mathrm{tet}} \rangle$")
        ax2.legend()
        ax2.grid(True, alpha=0.3)

        if has_openmm and t_o.size:
            _plot_tk(ax3, t_o, tk_o, "OpenMM CUDA + UMA", "b-")
        if reference_data is not None:
            t_r, _, _, tk_r = reference_data
            if t_r.size:
                _plot_tk(ax3, t_r, tk_r, "OpenMM Reference + UMA", "c-")
        if cuda_sequential_data is not None:
            t_s, _, _, tk_s = cuda_sequential_data
            if t_s.size:
                _plot_tk(ax3, t_s, tk_s, "OpenMM CUDA + UMA (sequential)", "m-.", lw=1.6)
        if has_ipi and t_i.size:
            _plot_tk(ax3, t_i, tk_i, "i-PI + LAMMPS UMA", "r--")
        if tip4p_data is not None:
            t_t, q6_t, qt_t, tk_t = tip4p_data
            _plot_tk(ax3, t_t, tk_t, "TIP4P/2005f RPMD", "g-.", lw=1.8)
        if not args.no_temperature_line:
            ax3.axhline(
                args.bath_temperature_k,
                color="k",
                ls=":",
                lw=1.0,
                alpha=0.6,
                label=f"Bath target ({args.bath_temperature_k:.0f} K)",
            )
        ax3.set_ylabel(r"$T_{\mathrm{K,centroid}}$ (K)")
        ax3.set_xlabel("Time (ps)")
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        tk_y_top = float(args.bath_temperature_k) * float(args.tk_ymax_factor)
        y_bottom, _ = ax3.get_ylim()
        ax3.set_ylim(bottom=max(0.0, float(y_bottom)), top=tk_y_top)
        if t_max is not None:
            ax3.set_xlim(0.0, float(t_max))

        fig.suptitle(_enrich_title_with_system(title_fmt, n_mol, n_beads_ann))
        # Tighter vertical packing between stacked axes (title spacing is handled above).
        fig.set_constrained_layout_pads(w_pad=0.02, h_pad=0.0, hspace=0.06)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(args.output, dpi=args.dpi, bbox_inches="tight", pad_inches=0.08)
    print(f"Saved {args.output}")

    if args.cv_plane is not None:
        _save_cv_plane_figure(
            args.cv_plane,
            has_openmm=has_openmm and bool(t_o.size),
            q6_o=q6_o,
            qt_o=qt_o,
            reference_data=reference_data,
            cuda_sequential_data=cuda_sequential_data,
            has_ipi=has_ipi and bool(t_i.size),
            q6_i=q6_i,
            qt_i=qt_i,
            tip4p_data=tip4p_data,
            title=args.title,
            n_molecules=n_mol,
            n_beads_ann=n_beads_ann,
            dpi=args.dpi,
        )
        print(f"Saved CV-plane figure {args.cv_plane}")


if __name__ == "__main__":
    main()
