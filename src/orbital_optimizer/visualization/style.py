"""Shared visual language of DAEDALUS figures: palette, ink and Matplotlib axes style.

One place for the colours used by the exported figures and by the
Verification Console, so a series has the same colour everywhere.

Rules applied by every figure (they are encoding rules, not decoration):
  - colour follows the entity, never its rank: an algorithm keeps its colour
    in every chart, whichever other algorithms are shown (``algorithm_color``);
  - categorical hues are assigned in the fixed order below and never cycled;
  - the null baseline (random search) is drawn in the de-emphasis grey;
  - one y-scale per plot; solid hairline grid; 2-point lines; text in ink
    colours, never in a series colour;
  - every figure has a table twin (the CSV/JSON written next to it).

Palette: eight categorical hues in an order that keeps adjacent pairs apart
under colour-vision deficiency. The slots used here were checked on
2026-11-10 with a palette validator (lightness band, chroma floor, adjacent
CVD and normal-vision separation): slots 1-5 pass on both surfaces; on the
light surface aqua, yellow and magenta have less than 3:1 contrast, which is
why identity is also carried by the legend, row labels and the table twin.
"""

from __future__ import annotations

from typing import Any

SURFACE: dict[str, str] = {"light": "#fcfcfb", "dark": "#1a1a19"}
INK: dict[str, dict[str, str]] = {
    "light": {"primary": "#0b0b0b", "secondary": "#52514e", "muted": "#898781",
              "grid": "#e1e0d9", "axis": "#c3c2b7"},
    "dark": {"primary": "#f0efec", "secondary": "#c3c2b7", "muted": "#898781",
             "grid": "#2c2c2a", "axis": "#383835"},
}
CATEGORICAL: dict[str, tuple[str, ...]] = {
    "light": ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"),
    "dark": ("#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#008300", "#9085e9", "#e66767"),
}
# Status colours mean good / warning / serious / critical; they are never used for a series.
STATUS: dict[str, str] = {"good": "#0ca30c", "warning": "#fab219", "serious": "#ec835a", "critical": "#d03b3b"}
BODY_FILL: dict[str, str] = {"light": "#e1e0d9", "dark": "#383835"}

# Fixed slot of each algorithm family; None is the de-emphasis grey of the null baseline.
ALGORITHM_SLOTS: dict[str, int | None] = {
    "DE": 0, "PSO": 1, "GA": 2, "Hybrid": 3, "NelderMead": 4, "RandomSearch": None,
}
ALGORITHM_LABELS: dict[str, str] = {
    "DE": "Differential Evolution", "PSO": "Particle swarm", "GA": "Genetic algorithm",
    "Hybrid": "Hybrid (global + Nelder-Mead)", "NelderMead": "Nelder-Mead", "RandomSearch": "Random search",
}


def algorithm_family(name: str) -> str:
    """Family of an algorithm name such as 'DE/rand/1/bin@1' or 'Hybrid(DE/... -> ...)@1'."""
    if not isinstance(name, str) or not name:
        raise ValueError("algorithm name must be a nonempty string")
    return name.split("(")[0].split("/")[0].split("@")[0]


def algorithm_color(name: str, mode: str = "light") -> str:
    """Colour of an algorithm by its family, the same in every figure.

    A family that is not in ALGORITHM_SLOTS gets slot 6; a figure showing
    several unknown families must assign them distinct slots itself.
    """
    if mode not in CATEGORICAL:
        raise ValueError("mode must be 'light' or 'dark'")
    slot = ALGORITHM_SLOTS.get(algorithm_family(name), 5)
    return INK[mode]["muted"] if slot is None else CATEGORICAL[mode][slot]


def algorithm_label(name: str) -> str:
    """Readable label of an algorithm family; the family itself when unknown."""
    family = algorithm_family(name)
    return ALGORITHM_LABELS.get(family, family)


def style_axes(axes: Any, *, grid: str = "y") -> None:
    """Apply the recessive chrome to Matplotlib axes: surface, hairline solid grid, ink colours.

    ``grid``: 'x', 'y', 'both' or 'none'.
    """
    ink = INK["light"]
    axes.set_facecolor(SURFACE["light"])
    for side in ("top", "right"):
        axes.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        axes.spines[side].set_color(ink["axis"])
        axes.spines[side].set_linewidth(0.8)
    axes.tick_params(colors=ink["muted"], labelsize=8, length=3, width=0.8)
    axes.xaxis.label.set_color(ink["secondary"])
    axes.yaxis.label.set_color(ink["secondary"])
    axes.xaxis.label.set_size(9)
    axes.yaxis.label.set_size(9)
    axes.title.set_color(ink["primary"])
    axes.title.set_size(10)
    axes.set_axisbelow(True)
    if grid != "none":
        axes.grid(True, axis=grid, color=ink["grid"], linewidth=0.8, linestyle="-")
    else:
        axes.grid(False)
