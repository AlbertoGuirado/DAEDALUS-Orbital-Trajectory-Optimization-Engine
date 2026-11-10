"""Figures of a multi-seed benchmark (RF-6A.06, RF-6A.08), exported with Matplotlib.

Every figure is drawn from the benchmark summary or its run rows, never from
hand-entered numbers, and has a table twin: ``runs.csv`` and
``convergence.csv`` hold the plotted values. Axes carry units; a curve that
aggregates seeds states its statistic (median and interquartile band).

Gaps are drawn on a logarithmic axis. A gap that is zero or negative has no
logarithm: it is drawn at ``LOG_FLOOR`` and the axis label says so.

Matplotlib is optional (extra ``plots``) and imported inside the functions.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np

from ..visualization.style import INK, SURFACE, algorithm_color, algorithm_label, style_axes

LOG_FLOOR: float = 1e-16


def _pyplot() -> Any:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def _floor(values: Any) -> np.ndarray:
    """Values for a log axis: non-positive entries drawn at LOG_FLOOR, +inf kept out of range."""
    array = np.array(values, dtype=float)
    return np.where(array > LOG_FLOOR, array, LOG_FLOOR)


def _grid(count: int, plt: Any, *, width: float = 4.6, height: float = 3.4) -> tuple[Any, list[Any]]:
    columns = min(count, 3)
    rows = int(np.ceil(count / columns))
    figure, axes = plt.subplots(rows, columns, figsize=(width * columns, height * rows), squeeze=False)
    figure.patch.set_facecolor(SURFACE["light"])
    flat = list(axes.flat)
    for extra in flat[count:]:
        extra.set_visible(False)
    return figure, flat[:count]


def _quantity(summary: Mapping[str, Any], scenario: str) -> tuple[str, bool]:
    """('gap' or 'objective', has_reference) for a scenario of the summary."""
    cells = summary["scenarios"][scenario]
    has_reference = all(cell["gap"] is not None for cell in cells.values())
    return ("gap" if has_reference else "objective"), has_reference


def _legend(figure: Any, handles: Sequence[Any], labels: Sequence[str]) -> None:
    legend = figure.legend(handles, labels, loc="upper center", ncol=min(len(labels), 5), frameon=False,
                           fontsize=8, bbox_to_anchor=(0.5, 1.0))
    for text in legend.get_texts():
        text.set_color(INK["light"]["secondary"])


def _caption(figure: Any, text: str) -> None:
    """One-line note under the panels: what the marks and the aggregation mean."""
    figure.text(0.5, 0.015, text, ha="center", va="bottom", fontsize=8, color=INK["light"]["muted"])


def save_convergence_figure(summary: Mapping[str, Any], references: Mapping[str, float | None],
                            units: Mapping[str, str], output: Path) -> Path:
    """Median best-so-far curve per algorithm with its interquartile band, one panel per scenario.

    y: best objective minus the scenario reference (the gap) when a
    reference exists, else the best objective; logarithmic. x: objective
    evaluations. ``references`` and ``units`` map scenario id to its
    reference value and objective unit.
    """
    plt = _pyplot()
    scenarios = list(summary["scenarios"])
    figure, axes = _grid(len(scenarios), plt)
    handles: dict[str, Any] = {}
    for axis, scenario in zip(axes, scenarios):
        reference = references.get(scenario)
        offset = 0.0 if reference is None else reference
        for cell in summary["scenarios"][scenario].values():
            band = cell["convergence"]
            if band is None:
                continue
            name = cell["algorithm_name"]
            color = algorithm_color(name)
            x = np.array(band["evaluations"], dtype=float)
            median = np.array(band["median"], dtype=float) - offset
            low = np.array(band["q25"], dtype=float) - offset
            high = np.array(band["q75"], dtype=float) - offset
            shown = np.isfinite(median)
            axis.fill_between(x[shown], _floor(low[shown]), _floor(np.where(np.isfinite(high), high, median)[shown]),
                              color=color, alpha=0.12, linewidth=0)
            line, = axis.plot(x[shown], _floor(median[shown]), color=color, linewidth=2,
                              solid_capstyle="round", solid_joinstyle="round")
            handles.setdefault(algorithm_label(name), line)
        axis.set_yscale("log")
        axis.set_xlabel("Objective evaluations")
        what = "gap to reference" if reference is not None else "best objective"
        axis.set_ylabel(f"Best-so-far {what} ({units[scenario]})")
        axis.set_title(scenario)
        style_axes(axis, grid="both")
    _legend(figure, list(handles.values()), list(handles))
    _caption(figure, "Convergence: median across seeds, interquartile band "
                     f"(values <= {LOG_FLOOR:g} drawn at {LOG_FLOOR:g})")
    figure.tight_layout(rect=(0, 0.07, 1, 0.94))
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=160, facecolor=figure.get_facecolor())
    plt.close(figure)
    return output


def save_final_distribution_figure(rows: Sequence[Mapping[str, Any]], units: Mapping[str, str],
                                   output: Path) -> Path:
    """Final gap (or objective) of every run, one row per algorithm, with median and interquartile range.

    Each small dot is one seed; the ringed dot is the median and the bar the
    interquartile range (inverted-CDF quantiles). Runs without a valid point
    are listed in the row label instead of being drawn. Logarithmic x axis.
    """
    plt = _pyplot()
    scenarios = list(dict.fromkeys(row["scenario"] for row in rows))
    figure, axes = _grid(len(scenarios), plt, height=3.0)
    ink = INK["light"]
    for axis, scenario in zip(axes, scenarios):
        selected = [row for row in rows if row["scenario"] == scenario]
        algorithms = list(dict.fromkeys(row["algorithm"] for row in selected))
        use_gap = all(row["gap"] is not None for row in selected)
        labels = []
        for position, algorithm in enumerate(algorithms):
            runs = [row for row in selected if row["algorithm"] == algorithm]
            name = runs[0]["algorithm_name"]
            values = np.array([row["gap"] if use_gap else row["f_best"] for row in runs], dtype=float)
            finite = np.sort(values[np.isfinite(values)])
            color = algorithm_color(name)
            y = len(algorithms) - 1 - position
            label = algorithm_label(name)
            if finite.size < values.size:
                label += f"\n({values.size - finite.size} of {values.size} without a valid point)"
            labels.append((y, label))
            if not finite.size:
                continue
            # Deterministic vertical jitter, so the figure is identical on every export.
            jitter = (np.arange(finite.size) % 7 - 3) * 0.045
            axis.scatter(_floor(finite), y + jitter, s=10, color=color, alpha=0.45, linewidths=0, zorder=2)
            quartiles = [finite[max(1, int(np.ceil(q * finite.size))) - 1] for q in (0.25, 0.5, 0.75)]
            axis.plot(_floor([quartiles[0], quartiles[2]]), [y, y], color=color, linewidth=2,
                      solid_capstyle="round", zorder=3)
            axis.scatter(_floor([quartiles[1]]), [y], s=46, color=color, edgecolors=SURFACE["light"],
                         linewidths=1.5, zorder=4)
        axis.set_yticks([y for y, _ in labels])
        axis.set_yticklabels([label for _, label in labels])
        axis.set_ylim(-0.6, len(algorithms) - 0.4)
        axis.set_xscale("log")
        what = "Final gap to reference" if use_gap else "Final objective"
        axis.set_xlabel(f"{what} ({units[scenario]})")
        axis.set_title(scenario)
        style_axes(axis, grid="x")
        for tick in axis.get_yticklabels():
            tick.set_color(ink["secondary"])
    _caption(figure, "One dot per seed; ringed dot = median, bar = interquartile range "
                     f"(values <= {LOG_FLOOR:g} drawn at {LOG_FLOOR:g})")
    figure.tight_layout(rect=(0, 0.08, 1, 1))
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=160, facecolor=figure.get_facecolor())
    plt.close(figure)
    return output


def save_success_figure(summary: Mapping[str, Any], output: Path) -> Path | None:
    """Success rate per algorithm with its 95 % Wilson interval, one panel per scenario with a criterion.

    Returns None, writing nothing, when no scenario declares a success threshold.
    """
    scenarios = [scenario for scenario, cells in summary["scenarios"].items()
                 if all(cell["success"] is not None for cell in cells.values())]
    if not scenarios:
        return None
    plt = _pyplot()
    figure, axes = _grid(len(scenarios), plt, height=2.8)
    ink = INK["light"]
    for axis, scenario in zip(axes, scenarios):
        cells = summary["scenarios"][scenario]
        count = len(cells)
        for position, cell in enumerate(cells.values()):
            y = count - 1 - position
            color = algorithm_color(cell["algorithm_name"])
            low, high = cell["success"]["wilson_95"]
            axis.plot([100 * low, 100 * high], [y, y], color=color, linewidth=2, solid_capstyle="round", zorder=3)
            axis.scatter([100 * cell["success"]["rate"]], [y], s=46, color=color, edgecolors=SURFACE["light"],
                         linewidths=1.5, zorder=4)
            axis.annotate(f"{cell['success']['successes']}/{cell['success']['runs']}",
                          (100 * high, y), xytext=(6, 0), textcoords="offset points", va="center",
                          fontsize=8, color=ink["secondary"])
        axis.set_yticks(range(count))
        axis.set_yticklabels([algorithm_label(cell["algorithm_name"]) for cell in reversed(list(cells.values()))])
        axis.set_ylim(-0.6, count - 0.4)
        axis.set_xlim(-2, 112)
        axis.set_xticks([0, 25, 50, 75, 100])
        axis.set_xlabel("Runs meeting the success criterion (%)")
        axis.set_title(scenario)
        style_axes(axis, grid="x")
        for tick in axis.get_yticklabels():
            tick.set_color(ink["secondary"])
    _caption(figure, "Dot = observed rate, bar = 95 % Wilson interval, label = successes / seeds")
    figure.tight_layout(rect=(0, 0.09, 1, 1))
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=160, facecolor=figure.get_facecolor())
    plt.close(figure)
    return output


def save_cost_figure(summary: Mapping[str, Any], units: Mapping[str, str], output: Path) -> Path:
    """Median final gap against median wall time per run, one labelled point per algorithm (RF-6A.06).

    Bars span the interquartile ranges. Wall time depends on the machine and
    on concurrent load; it is shown as context and is never a stopping
    criterion. Both axes are logarithmic.
    """
    plt = _pyplot()
    scenarios = list(summary["scenarios"])
    figure, axes = _grid(len(scenarios), plt)
    ink = INK["light"]
    for axis, scenario in zip(axes, scenarios):
        quantity, _ = _quantity(summary, scenario)
        for cell in summary["scenarios"][scenario].values():
            stats, time = cell[quantity], cell["wall_time_s"]
            if not np.isfinite(stats["median"]):
                continue
            color = algorithm_color(cell["algorithm_name"])
            x, y = time["median"], float(_floor(stats["median"]))
            axis.plot([time["q1"], time["q3"]], [y, y], color=color, linewidth=2, solid_capstyle="round", zorder=3)
            upper = stats["q3"] if np.isfinite(stats["q3"]) else stats["median"]
            axis.plot([x, x], _floor([stats["q1"], upper]), color=color, linewidth=2, solid_capstyle="round", zorder=3)
            axis.scatter([x], [y], s=46, color=color, edgecolors=SURFACE["light"], linewidths=1.5, zorder=4)
            axis.annotate(algorithm_label(cell["algorithm_name"]), (x, y), xytext=(7, 5),
                          textcoords="offset points", fontsize=8, color=ink["secondary"])
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_xlabel("Wall time per run (s)")
        what = "Final gap to reference" if quantity == "gap" else "Final objective"
        axis.set_ylabel(f"{what} ({units[scenario]})")
        axis.set_title(scenario)
        style_axes(axis, grid="both")
    _caption(figure, "Dot = medians across seeds, bars = interquartile ranges; equal evaluation budgets "
                     f"(values <= {LOG_FLOOR:g} drawn at {LOG_FLOOR:g})")
    figure.tight_layout(rect=(0, 0.07, 1, 1))
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=160, facecolor=figure.get_facecolor())
    plt.close(figure)
    return output
