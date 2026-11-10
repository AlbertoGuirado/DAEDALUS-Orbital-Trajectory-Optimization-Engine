"""Figure of an impulsive transfer in its orbital plane (RF-6A.06), exported with Matplotlib.

The figure is drawn from a propagated ``ImpulsiveTrajectory``; it does not
compute physics of its own beyond projecting positions onto the plane of the
initial orbit and converting metres to kilometres at this presentation
boundary.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from ..simulation.impulsive import ImpulsiveTrajectory
from .style import BODY_FILL, CATEGORICAL, INK, SURFACE, style_axes


def orbital_plane_basis(state: NDArray[np.float64]) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Unit vectors (e1, e2) spanning the plane of motion of an SI state (6,).

    e1 is the radial direction of the state and e2 = h_hat x e1 the transverse
    direction (prograde), with h = r x v. A radial state has no plane and
    raises ValueError.
    """
    position, velocity = np.asarray(state[:3], dtype=float), np.asarray(state[3:], dtype=float)
    momentum = np.cross(position, velocity)
    norm = float(np.linalg.norm(momentum))
    if norm == 0 or not np.isfinite(norm):
        raise ValueError("the state is radial: its plane of motion is undefined")
    first = position / np.linalg.norm(position)
    return first, np.cross(momentum / norm, first)


def plane_coordinates(positions: NDArray[np.float64], basis: tuple[NDArray[np.float64], NDArray[np.float64]]
                      ) -> NDArray[np.float64]:
    """Project positions (M, 3), in metres, onto a plane basis; returns (M, 2) in metres."""
    return np.column_stack((positions @ basis[0], positions @ basis[1]))


def save_transfer_figure(trajectory: ImpulsiveTrajectory, initial_radius: float, final_radius: float,
                         output: Path, *, body_radius: float | None = None, title: str | None = None) -> Path:
    """Save the initial orbit, the transfer arcs, the target orbit and the burn points; axes in km.

    ``trajectory`` should be sampled densely (``samples_per_arc``). Radii in
    metres; ``body_radius`` draws the central body to scale when given. The
    positions are projected onto the plane of the first sampled state, which
    is the plane of the whole transfer for coplanar burns.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ink, colors = INK["light"], CATEGORICAL["light"]
    states = trajectory.states()
    basis = orbital_plane_basis(trajectory.arcs[0].states[0])
    angle = np.linspace(0.0, 2 * np.pi, 721)
    figure, axis = plt.subplots(figsize=(6.4, 6.4))
    figure.patch.set_facecolor(SURFACE["light"])
    if body_radius is not None:
        axis.fill(body_radius / 1000 * np.cos(angle), body_radius / 1000 * np.sin(angle),
                  color=BODY_FILL["light"], linewidth=0, zorder=1)
    axis.plot(initial_radius / 1000 * np.cos(angle), initial_radius / 1000 * np.sin(angle),
              color=colors[0], linewidth=1.2, label="Initial orbit", zorder=2)
    axis.plot(final_radius / 1000 * np.cos(angle), final_radius / 1000 * np.sin(angle),
              color=colors[2], linewidth=1.2, label="Target orbit", zorder=2)
    for index, arc in enumerate(trajectory.arcs):
        path = plane_coordinates(arc.states[:, :3], basis) / 1000
        axis.plot(path[:, 0], path[:, 1], color=colors[1], linewidth=2, solid_capstyle="round",
                  label="Transfer trajectory" if index == 0 else None, zorder=3)
    # Position of each burn: the sample of the trajectory at the burn time.
    times = trajectory.times()
    for number, (time, vector) in enumerate(zip(trajectory.burn_times, trajectory.burn_vectors), start=1):
        point = plane_coordinates(states[int(np.argmin(np.abs(times - time)))][None, :3], basis)[0] / 1000
        axis.scatter([point[0]], [point[1]], s=52, color=ink["primary"], edgecolors=SURFACE["light"],
                     linewidths=1.5, zorder=5)
        axis.annotate(f"Δv{number} = {np.linalg.norm(vector):.1f} m/s", point, xytext=(8, 8),
                      textcoords="offset points", fontsize=8, color=ink["secondary"], zorder=6)
    extent = 1.08 * max(float(np.max(np.linalg.norm(states[:, :3], axis=1))), initial_radius, final_radius) / 1000
    axis.set_xlim(-extent, extent)
    axis.set_ylim(-extent, extent)
    axis.set_aspect("equal")
    axis.set_xlabel("Position along the departure radius (km)")
    axis.set_ylabel("Position along the departure velocity (km)")
    if title:
        axis.set_title(title)
    style_axes(axis, grid="both")
    legend = axis.legend(loc="upper right", frameon=False, fontsize=8)
    for text in legend.get_texts():
        text.set_color(ink["secondary"])
    figure.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=160, facecolor=figure.get_facecolor())
    plt.close(figure)
    return output
